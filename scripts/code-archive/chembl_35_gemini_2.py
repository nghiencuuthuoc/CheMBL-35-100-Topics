import os
import pandas as pd
from google import genai
from dotenv import dotenv_values
from markdown import markdown
import subprocess

# === CONFIGURATION ===
# excel_path = "../prompt/chembl_35_topics_bilingual_1_100_prompt.xlsx"
excel_path = "../prompt/prompt_test.xlsx"
output_folder = "../chembl-35-100-topics-ai"
template_docx = "../template/template.docx"
env_path = "../env/.env"
model_name = "gemini-2.0-flash"

# === LOAD API KEYS FROM .env FILE ===
env_vars = dotenv_values(env_path)
api_keys = [v for k, v in env_vars.items() if k == "GEMINI_API_KEY"]
if not api_keys:
    raise ValueError("❌ No GEMINI_API_KEY found in the .env file.")

# === READ PROMPTS FROM EXCEL ===
os.makedirs(output_folder, exist_ok=True)
df = pd.read_excel(excel_path, engine="openpyxl")

# === TRY GEMINI API WITH FALLBACK ON MULTIPLE KEYS ===
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
                print(f"⚠️  503 error with key {key[:20]}... switching to the next key.")
                continue
            else:
                raise e
    raise RuntimeError("❌ All API keys failed. Please check your .env file.")

# === LOOP THROUGH PROMPTS AND GENERATE OUTPUT ===
for idx, row in df.iterrows():
    prompt = "\n".join(str(cell) for cell in row if pd.notna(cell))
    try:
        output_text = try_generate(prompt, api_keys)

        base_filename = f"topic_{61 + idx:03d}"
        md_path = os.path.join(output_folder, f"{base_filename}.md")
        html_path = os.path.join(output_folder, f"{base_filename}.html")
        docx_path = os.path.join(output_folder, f"{base_filename}.docx")

        # Save markdown file
        with open(md_path, "w", encoding="utf-8") as f_md:
            f_md.write(output_text)

        # Save HTML file
        with open(html_path, "w", encoding="utf-8") as f_html:
            f_html.write(markdown(output_text))

        # Convert markdown to DOCX using pandoc and custom template
        subprocess.run([
            "pandoc", md_path,
            "-o", docx_path,
            "--reference-doc", template_docx
        ])

        print(f"✅ Generated: {base_filename}")

    except Exception as e:
        print(f"❌ Error at row {idx + 1}: {e}")
