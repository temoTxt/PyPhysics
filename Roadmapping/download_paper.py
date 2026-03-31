import arxiv
import os
import re

def sanitize_filename(title):
    """Removes characters that are illegal in file names across OSs."""
    return re.sub(r'[\\/*?:"<>|]', "", title).strip()

def download_tepper_gill_papers(download_path="./Tepper_Gill_Papers"):
    # 1. Initialize the Client (The modern way to handle requests)
    client = arxiv.Client(
        page_size=100,
        delay_seconds=3,
        num_retries=3
    )

    # 2. Define the search
    search = arxiv.Search(
        query='au:"Tepper L. Gill"',
        sort_by=arxiv.SortCriterion.SubmittedDate
    )

    if not os.path.exists(download_path):
        os.makedirs(download_path)

    print(f"Starting downloads for Tepper L. Gill via arxiv.Client...\n")

    # 3. Use client.results(search) instead of search.results()
    count = 0
    for result in client.results(search):
        clean_title = sanitize_filename(result.title)
        # Adding .pdf extension to the filename
        filename = f"{clean_title}.pdf"
        
        print(f"Downloading [{count+1}]: {result.title[:50]}...")
        
        try:
            # This downloads the latest version by default
            result.download_pdf(dirpath=download_path, filename=filename)
            count += 1
        except Exception as e:
            print(f"Error downloading {result.title}: {e}")

    print(f"\nFinished! {count} papers saved to '{download_path}'.")

if __name__ == "__main__":
    download_tepper_gill_papers()