"""Inspect the public DDSE interface in an isolated browser; no user profile."""
from pathlib import Path
import sys, json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/vendor'))
from playwright.sync_api import sync_playwright

def main():
    dest=ROOT/'data/raw/unified_sources';dest.mkdir(exist_ok=True)
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless=True)
        context=browser.new_context(accept_downloads=True)
        page=context.new_page()
        events=[]
        page.on('response',lambda r: events.append({'url':r.url.split('?')[0],'status':r.status}) if 'digbat' in r.url or 'ddse' in r.url else None)
        try:
            page.goto('https://www.digbat.org/',wait_until='domcontentloaded',timeout=60000)
            page.wait_for_timeout(18000)
            body=page.locator('body').inner_text()
            (dest/'ddse_browser_body.txt').write_text(body)
            (dest/'ddse_browser.html').write_text(page.content())
            page.screenshot(path=str(dest/'ddse_browser.png'),full_page=True)
            (dest/'ddse_browser_events.json').write_text(json.dumps(events,indent=2))
            print(body[:7000])
            print('buttons:',page.get_by_role('button').all_text_contents())
            print('links:',page.locator('a').evaluate_all('(els)=>els.map(a=>({text:a.innerText,url:a.href}))'))
        except Exception as e: print(type(e).__name__,str(e)[:500])
        browser.close()
if __name__=='__main__':main()
