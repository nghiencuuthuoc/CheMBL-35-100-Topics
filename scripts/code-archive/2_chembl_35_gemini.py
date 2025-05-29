import os
import time
import pandas as pd
from google import genai
from dotenv import load_dotenv

# Load API key from .env file
load_dotenv("../env/.env")

# Set up Gemini client and model
client = genai.Client(api_key=os.getenv("GEMINI_API_KEY"))
model = client.GenerativeModel("gemini-1.5-flash")

# Load Excel file
file_path = "../data/chembl_35_topics_bilingual_61_100.xlsx"
df = pd.read_excel(file_path)

# Define output folder
output_folder = "../chembl-35-100-topics-ai"
os.makedirs(output_folder, exist_ok=True)

# Loop through each row
for index, row in df.iterrows():
    try:
        # Combine all non-null cells into a single prompt
        prompt = "\n".join([str(cell).strip() for cell in row if pd.notna(cell)])

        retries = 3
        for attempt in range(retries):
            try:
                response = model.generate_content(prompt)
                output = response.text.strip()

                # Save to .md
                md_file = os.path.join(output_folder, f"topic_{index+61:03d}.md")
                with open(md_file, "w", encoding="utf-8") as f_md:
                    f_md.write(output)

                # Save to .html
                html_file = os.path.join(output_folder, f"topic_{index+61:03d}.html")
                with open(html_file, "w", encoding="utf-8") as f_html:
                    f_html.write(f"<html><body><pre>{output}</pre></body></html>")

                print(f"✅ Generated: topic_{index+61:03d}")
                break
            except Exception as e:
                if "503" in str(e) and attempt < retries - 1:
                    print(f"⚠️ Retry {attempt+1}/3 after 5s due to 503 overload.")
                    time.sleep(5)
                else:
                    raise e
    except KeyboardInterrupt:
        print("❌ Interrupted by user.")
        break
    except Exception as e:
        print(f"❌ Error at row {index+61}: {e}")
