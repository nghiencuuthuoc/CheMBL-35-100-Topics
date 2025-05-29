# chembl_35_topics_bilingual_1_100_prompt

import os
import time
import pandas as pd
import markdown
from weasyprint import HTML
from google import genai
from dotenv import load_dotenv
import pypandoc  # Thêm để xuất .docx

# Load API key từ .env
load_dotenv("../env/.env")
client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))

# Input/output
file_path = "../prompt/chembl_35_topics_100_prompt.xlsx"
output_folder = "../chembl-35-100-topics-ai-en"
os.makedirs(output_folder, exist_ok=True)

# Đọc file Excel, bỏ dòng tiêu đề đầu
df = pd.read_excel(file_path, header=1)

# Lặp qua từng dòng
for index, row in df.iterrows():
    try:
        prompt = "\n".join([str(cell).strip() for cell in row if pd.notna(cell)])
        retries = 3
        for attempt in range(retries):
            try:
                response = client.models.generate_content(
                    model="gemini-2.0-flash",
                    contents=prompt
                )
                output = response.text.strip()

                # Tạo nội dung markdown
                combined_md = f"{prompt}\n\n---\n\n{output}"
                base_name = f"{index+1:03d}_chembl_topic"
                md_path = os.path.join(output_folder, f"{base_name}.md")
                html_path = os.path.join(output_folder, f"{base_name}.html")
                pdf_path = os.path.join(output_folder, f"{base_name}.pdf")
                docx_path = os.path.join(output_folder, f"{base_name}.docx")

                # Lưu file .md
                with open(md_path, "w", encoding="utf-8") as f_md:
                    f_md.write(combined_md)

                # Chuyển Markdown sang HTML
                html_body = markdown.markdown(combined_md, output_format="html5")
                full_html = f"<html><body><pre>{html_body}</pre></body></html>"

                # Lưu .html
                with open(html_path, "w", encoding="utf-8") as f_html:
                    f_html.write(full_html)

                # Lưu .pdf qua WeasyPrint
                HTML(string=full_html).write_pdf(pdf_path)

                # Lưu .docx qua pypandoc
                pypandoc.convert_file(md_path, 'docx', outputfile=docx_path)

                print(f"✅ Generated: {base_name}")
                break
            except Exception as e:
                if "503" in str(e) and attempt < retries - 1:
                    print(f"⚠️ Retry {attempt+1}/3 after 5s due to 503.")
                    time.sleep(5)
                else:
                    raise e
    except KeyboardInterrupt:
        print("⛔ Stopped by user.")
        break
    except Exception as e:
        print(f"❌ Error at row {index+1}: {e}")
