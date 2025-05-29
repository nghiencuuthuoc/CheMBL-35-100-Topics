import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk
import pandas as pd
from markdown import markdown
from dotenv import dotenv_values
from google import genai

# Header markdown
header_md = (
    "# PharmApp Suite\n"
    "## 🧠 AI for Drug Discovery and Development 🧪\n"
    "| Copyright 2025 | Nghiên Cứu Thuốc | www.nghiencuuthuoc.com | Zalo: +84888999311 |\n"
)

# Load API keys
def load_api_keys(env_path):
    env_vars = dotenv_values(env_path)
    keys = [v for k, v in env_vars.items() if k == "GEMINI_API_KEY"]
    if not keys:
        raise ValueError("❌ No GEMINI_API_KEY found in .env")
    return keys

# Try generating content using fallback keys
def try_generate(prompt, keys, model_name):
    for key in keys:
        try:
            client = genai.Client(api_key=key)
            response = client.models.generate_content(model=model_name, contents=prompt)
            return response.text
        except Exception as e:
            if "503" in str(e):
                print(f"⚠️  503 error with key {key[:20]}... switching.")
                continue
            else:
                raise e
    raise RuntimeError("❌ All API keys failed.")

# Main processing logic
def run_generation():
    try:
        excel_path = prompt_file.get()
        output_folder = output_dir.get()
        env_path = env_file.get()
        template_docx = template_path.get()
        soffice_path = soffice_exec.get()
        model_name = "gemini-2.0-flash"

        os.makedirs(output_folder, exist_ok=True)
        keys = load_api_keys(env_path)
        df = pd.read_excel(excel_path, engine="openpyxl")

        for idx, row in df.iterrows():
            try:
                title = str(row.iloc[4]) if pd.notna(row.iloc[4]) else "Untitled"
                title = '🧩 Topic: ' + title
                prompt_parts = [str(cell) for i, cell in enumerate(row) if pd.notna(cell) and i != 4]
                prompt = "\n".join(prompt_parts)

                output_text = try_generate(prompt, keys, model_name)
                base_filename = str(row.iloc[2]) if pd.notna(row.iloc[2]) else f"unknown_{idx+1}"
                md_path = os.path.join(output_folder, f"{base_filename}.md")
                html_path = os.path.join(output_folder, f"{base_filename}.html")
                docx_path = os.path.join(output_folder, f"{base_filename}.docx")
                pdf_path = os.path.join(output_folder, f"{base_filename}.pdf")

                with open(md_path, "w", encoding="utf-8") as f_md:
                    f_md.write(f"{header_md}\n{title}\n---\n{output_text}")

                html_body = markdown(
                    f"{header_md}\n\n# {title}\n\n---\n\n{output_text}",
                    output_format="html5"
                )
                html_full = f"""<!DOCTYPE html><html><head><meta charset="UTF-8"><title>{base_filename}</title></head>
                <body style="font-family:Georgia; max-width:800px; margin:auto; padding:2em;">{html_body}</body></html>"""

                with open(html_path, "w", encoding="utf-8") as f_html:
                    f_html.write(html_full)

                subprocess.run([
                    "pandoc", md_path,
                    "-o", docx_path,
                    "--reference-doc", template_docx
                ], check=True)

                subprocess.run([
                    soffice_path, "--headless", "--convert-to", "pdf", os.path.basename(docx_path),
                    "--outdir", "."
                ], cwd=os.path.dirname(docx_path), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                status_var.set(f"✅ Done: {base_filename}")
                root.update()

            except Exception as e:
                print(f"❌ Error at row {idx + 1}: {e}")
                status_var.set(f"❌ Error at row {idx + 1}")
                root.update()
        messagebox.showinfo("Hoàn tất", "Tất cả các file đã xử lý xong.")
    except Exception as e:
        messagebox.showerror("Lỗi", str(e))

# GUI setup
root = tk.Tk()
root.title("🧠 ChEMBL 35 Generator GUI")
root.geometry("850x600")
root.configure(bg="white")

def browse_file(var, filetypes):
    path = filedialog.askopenfilename(filetypes=filetypes)
    if path:
        var.set(path)

def browse_folder(var):
    path = filedialog.askdirectory()
    if path:
        var.set(path)

def create_entry_row(label, var, row, browse_type=None, filetypes=None):
    tk.Label(frame, text=label, bg="white").grid(row=row, column=0, sticky="w")
    tk.Entry(frame, textvariable=var, width=70).grid(row=row, column=1)
    if browse_type == "file":
        tk.Button(frame, text="Browse", font=("Arial", 10), command=lambda: browse_file(var, filetypes)).grid(row=row, column=2)
    elif browse_type == "folder":
        tk.Button(frame, text="Browse", font=("Arial", 10), command=lambda: browse_folder(var)).grid(row=row, column=2)

# Logo
logo_path = "images/nct_logo.png"
if os.path.exists(logo_path):
    logo_img = Image.open(logo_path).resize((120, 120))
    logo = ImageTk.PhotoImage(logo_img)
    tk.Label(root, image=logo, bg="white").pack(pady=10)

tk.Label(root, text="Công cụ sinh tài liệu nghiên cứu thuốc từ Excel Prompt + Gemini", font=("Arial", 12), bg="white").pack()

frame = tk.Frame(root, bg="white", padx=20, pady=20)
frame.pack()

# Variables
prompt_file = tk.StringVar()
output_dir = tk.StringVar()
env_file = tk.StringVar()
template_path = tk.StringVar()
soffice_exec = tk.StringVar()

# Entry rows
create_entry_row("📄 Excel Prompt:", prompt_file, 0, "file", [("Excel files", "*.xlsx")])
create_entry_row("📂 Output Folder:", output_dir, 1, "folder")
create_entry_row("🔑 .env File:", env_file, 2, "file", [("Env files", "*.env")])
create_entry_row("📄 Template DOCX:", template_path, 3, "file", [("DOCX files", "*.docx")])
create_entry_row("📎 LibreOffice soffice.exe:", soffice_exec, 4, "file", [("Executable", "*.exe")])

# Run button
tk.Button(root, text="🚀 Generate Documents", font=("Arial", 14, "bold"), bg="#4CAF50", fg="white",
          padx=20, pady=10, command=run_generation).pack(pady=20)

status_var = tk.StringVar()
tk.Label(root, textvariable=status_var, bg="white", fg="blue", font=("Arial", 10)).pack()

# Footer
tk.Label(root, text="| Copyright 2025 | 🧠 Nghiên Cứu Thuốc | PharmApp |\n"
                    "| www.nghiencuuthuoc.com | Zalo: +84888999311 |",
         font=("Arial", 10), bg="white", fg="gray").pack(side="bottom", pady=10)

root.mainloop()
