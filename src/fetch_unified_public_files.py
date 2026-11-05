"""Download public source documents using an isolated Chrome request context."""
from pathlib import Path
import sys, json, hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/vendor'))
from playwright.sync_api import sync_playwright
def main():
    dest=ROOT/'data/raw/unified_sources'
    with sync_playwright() as p:
        b=p.chromium.launch(executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless=True)
        c=b.new_context()
        page=c.new_page()
        try:
            page.goto('https://authors.library.caltech.edu/records/64j4e-7fx68',wait_until='domcontentloaded',timeout=30000)
            print('Caltech page:',page.title(),page.locator('body').inner_text()[:180])
            if page.locator('a[href*="d2ee03499a1.pdf"]').count():
                with page.expect_download(timeout=30000) as info:
                    page.locator('a[href*="d2ee03499a1.pdf"]').first.click()
                info.value.save_as(str(dest/'laskowski2023_si.pdf'))
                print('Saved Laskowski public download')
        except Exception as e:print('browser',type(e).__name__,str(e).split('https:')[0][:80])
        jobs=[('https://authors.library.caltech.edu/records/64j4e-7fx68/files/d2ee03499a1.pdf?download=1','laskowski2023_si.pdf'),
              ('https://pcwww.liverpool.ac.uk/~msd30/lmds/LiIonDatabase.csv','liion_original.csv')]
        for url,name in jobs:
            try:
                r=c.request.get(url,timeout=60000)
                raw=r.body()
                valid=r.ok and (raw.startswith(b'%PDF') if name.endswith('.pdf') else b'composition' in raw[:400].lower())
                print(name,r.status,len(raw),'valid',valid)
                if valid:(dest/name).write_bytes(raw)
            except Exception as e:print(name,type(e).__name__,str(e).split('https:')[0][:100])
        b.close()
if __name__=='__main__':main()
