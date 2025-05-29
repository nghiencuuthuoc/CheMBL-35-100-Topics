# ok
import os
import pandas as pd
from google import genai
from dotenv import dotenv_values
from markdown import markdown
import subprocess

# === CONFIGURATION ===
excel_path = "../prompt/chembl_35_topics_bilingual_1_100_prompt.xlsx"
# excel_path = "../prompt/prompt_test.xlsx"
output_folder = "../chembl-35-100-topics-ai"
template_docx = "../template/template.docx"
env_path = "../env/.env"
model_name = "gemini-2.0-flash"

# Path to LibreOffice Portable soffice.exe
soffice_path = r"E:\PharmAppDev\PharmAppSuite-v2025.02\Apps\LibreOffice\program\soffice.exe"

# === LOAD API KEYS FROM .env ===
env_vars = dotenv_values(env_path)
api_keys = [v for k, v in env_vars.items() if k == "GEMINI_API_KEY"]
if not api_keys:
    raise ValueError("❌ No GEMINI_API_KEY found in .env")

# === PREPARE OUTPUT FOLDER ===
os.makedirs(output_folder, exist_ok=True)
df = pd.read_excel(excel_path, engine="openpyxl")

# === FUNCTION: Generate content with fallback API keys ===
def try_generate(prompt, keys):
    for key in keys:
        try:
            client = genai.Client(api_key=key)
            response = client.models.generate_content(
                model=model_name,
                contents=prompt
            )
            return response.text
        except Exception as e:
            if "503" in str(e):
                print(f"⚠️  503 error with key {key[:20]}... switching to next key.")
                continue
            else:
                raise e
    raise RuntimeError("❌ All API keys failed.")

# === MAIN LOOP: Generate all output files ===
for idx, row in df.iterrows():
    prompt = "\n".join(str(cell) for cell in row if pd.notna(cell))
    try:
        output_text = try_generate(prompt, api_keys)

        base_filename = f"topic_{61 + idx:03d}"
        md_path = os.path.join(output_folder, f"{base_filename}.md")
        html_path = os.path.join(output_folder, f"{base_filename}.html")
        docx_path = os.path.join(output_folder, f"{base_filename}.docx")
        pdf_path = os.path.join(output_folder, f"{base_filename}.pdf")

        # Save Markdown
        with open(md_path, "w", encoding="utf-8") as f_md:
            f_md.write(output_text)

        # Generate styled HTML
        html_body = markdown(output_text)
        html_full = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>{base_filename}</title>
            <style>
                body {{
                    font-family: "Georgia", serif;
                    max-width: 800px;
                    margin: 40px auto;
                    padding: 20px;
                    line-height: 1.6;
                    background-color: #ffffff;
                    color: #333;
                }}
                h1, h2, h3 {{
                    color: #1a1a1a;
                }}
                code {{
                    background-color: #f5f5f5;
                    padding: 2px 4px;
                    border-radius: 4px;
                }}
                pre {{
                    background-color: #f5f5f5;
                    padding: 10px;
                    overflow-x: auto;
                    border-radius: 6px;
                }}
            </style>
        </head>
        <body>
        {html_body}
        </body>
        </html>
        """
        with open(html_path, "w", encoding="utf-8") as f_html:
            f_html.write(html_full)

        # Convert to DOCX using Pandoc
        subprocess.run([
            "pandoc", md_path,
            "-o", docx_path,
            "--reference-doc", template_docx
        ], check=True)

        # Convert DOCX to PDF using LibreOffice (run from .docx folder)
        docx_dir = os.path.dirname(docx_path)
        docx_filename = os.path.basename(docx_path)

        result = subprocess.run([
            soffice_path, "--headless", "--convert-to", "pdf", docx_filename, "--outdir", "."
        ], cwd=docx_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        print("📤 LibreOffice stdout:")
        print(result.stdout)
        print("⚠️ LibreOffice stderr:")
        print(result.stderr)

        # Check PDF existence
        if not os.path.exists(pdf_path):
            print(f"❌ PDF not found: {pdf_path}")
            found = [f for f in os.listdir(output_folder) if f.endswith(".pdf")]
            print("📁 PDFs in folder:", found)
        else:
            print(f"✅ PDF created: {pdf_path}")

        print(f"✅ Generated: {base_filename}")

    except Exception as e:
        print(f"❌ Error at row {idx + 1}: {e}")
