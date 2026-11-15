// Headless Mol* renders for step 17. Reads scenes.json (written by
// 17_mechanism_cards.py), renders every fold twice (full, zoom) in Chromium
// with page.js, and writes transparent PNGs plus views.json, which records
// where the zoom window and the labels landed so the cards can draw on top.
//
// Usage: node render.mjs <scenes.json> <outdir>

import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"
import { chromium } from "playwright-core"

const HERE = path.dirname(fileURLToPath(import.meta.url))
const [scenesPath, outdir] = process.argv.slice(2)
const scenes = JSON.parse(fs.readFileSync(scenesPath, "utf8"))

// the Chromium Playwright already downloaded on this machine, or Chrome
const executablePath = (() => {
  const cache = path.join(process.env.HOME, "Library/Caches/ms-playwright")
  const found = fs.existsSync(cache)
    ? fs.readdirSync(cache).filter((d) => /^chromium-\d+$/.test(d)).sort().reverse()
    : []
  for (const d of found) {
    for (const app of ["Chromium.app/Contents/MacOS/Chromium",
      "Google Chrome for Testing.app/Contents/MacOS/Google Chrome for Testing"]) {
      const p = path.join(cache, d, "chrome-mac-arm64", app)
      if (fs.existsSync(p)) return p
    }
  }
  return "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
})()

const SIZE = 900          // css px; screenshots at 4x, 3600 px
const PAGE = `<!doctype html><html><head><meta charset="utf-8">
<style>
  html, body { margin: 0; background: transparent; }
  #frame { position: relative; width: ${SIZE}px; height: ${SIZE}px; }
  #app { position: absolute; inset: 0; }
  #labels { position: absolute; inset: 0; pointer-events: none; }
  .msp-layout-top, .msp-layout-bottom, .msp-layout-left, .msp-layout-right,
  .msp-layout-controls, .msp-viewport-controls, .msp-highlight-toast-wrapper,
  .msp-viewport-top-left-controls { display: none !important; }
  .msp-plugin, .msp-plugin-content, .msp-layout-standard, .msp-layout-expanded,
  .msp-viewport { background: transparent !important; position: absolute; inset: 0; }
  #app * { border: 0 !important; outline: 0 !important; box-shadow: none !important; }
  .label { position: absolute; transform: translate(-50%, -50%); white-space: nowrap;
    font: 600 28px/1 "Helvetica Neue", Arial, sans-serif; color: #1a1a1a;
    -webkit-text-stroke: 5px rgba(255,255,255,0.92); paint-order: stroke fill; }
  .label.secondary { font-weight: 400; font-size: 24px; color: #7a5a10; }
  .label.primary { color: #a3245f; }
  .label.substrate { color: #136b5d; font-style: italic; font-weight: 600; font-size: 24px; }
  .label.cofactor { color: #3c444c; font-weight: 600; font-size: 22px; }
</style></head><body><div id="frame"><div id="app"></div><div id="labels"></div></div></body></html>`

const browser = await chromium.launch({
  executablePath,
  args: ["--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
})
const page = await browser.newPage({ viewport: { width: SIZE, height: SIZE }, deviceScaleFactor: 4 })
page.on("pageerror", (e) => console.error("page error:", e.message))
await page.setContent(PAGE)
await page.addStyleTag({ path: path.join(HERE, "vendor/molstar.css") })
await page.addScriptTag({ path: path.join(HERE, "vendor/molstar.js") })
await page.addScriptTag({ path: path.join(HERE, "page.js") })

// a rerun of some scenes keeps the views of the others
const viewsPath = path.join(outdir, "views.json")
const views = fs.existsSync(viewsPath) ? JSON.parse(fs.readFileSync(viewsPath, "utf8")) : {}
for (const s of scenes) {
  views[s.slug] = {}
  for (const v of ["full", "zoom"]) {
    const zoom = v === "zoom"
    const spec = { ...s, pdb: fs.readFileSync(s.pdb_path, "utf8") }
    const res = await page.evaluate(
      ([spec, view]) => window.renderScene(spec, view),
      [spec, { zoom, labels: zoom || s.full_label_positions === true,
        hideLabels: !zoom || s.zoom_labels === false, size: SIZE, zoomRadius: s.zoom_radius, tighten: 1.08 }])
    fs.mkdirSync(path.join(outdir, v), { recursive: true })
    await page.locator("#frame").screenshot({ path: path.join(outdir, v, `${s.slug}.png`), omitBackground: true })
    views[s.slug][v] = res
  }
  console.log(`  ${s.slug}`)
}
fs.writeFileSync(viewsPath, JSON.stringify(views, null, 1))
await browser.close()
