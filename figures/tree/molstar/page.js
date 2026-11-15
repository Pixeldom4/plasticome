/*
 * Browser side of the Mol* renders (step 17). Loaded into a blank page by
 * render.mjs after vendor/molstar.js, the prebuilt viewer the iGEM wiki
 * serves. The drawing follows the wiki's StructureViewer
 * (iGEM_wiki/src/components/mdx/interactive/StructureViewer.js):
 *
 *   cartoon           Mol*'s own cartoon from its secondary structure
 *                     assignment, in the fold's colour; a second protein
 *                     chain lighter
 *   catalytic         pink sticks
 *   secondary         amber, thinner sticks, so they read as support rather
 *                     than the site itself
 *   substrate         chain L, teal balls at 70% of van der Waals size: a
 *                     bound molecule, not a residue (sticks for residues,
 *                     balls for the substrate)
 *   cofactors         heme as dark slate sticks, metal ions as silver balls
 *
 * with Mol*'s default lighting and occlusion. Added for print: an
 * orthographic camera set from a direction the Python side computes (the
 * active site faces the reader, the long axis of the fold lies across), and
 * in the zoom, labels for residues, the plastic and the cofactors, drawn as
 * HTML over the canvas at their projected positions so they stay sharp.
 */

const COLORS = {
  catalytic: 0xe8559a,
  secondary: 0xe0a526,
  substrate: 0x1f9e89,
  heme: 0x4a5560,
  metal: 0xb4b9bf,
}

const hex = (v) => (typeof v === "string" ? parseInt(v.replace("#", ""), 16) : v)
const wait = (ms) => new Promise((r) => setTimeout(r, ms))

let viewer = null

async function boot(background) {
  viewer = await molstar.Viewer.create(document.getElementById("app"), {
    extensions: [],
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    viewportShowControls: false,
    viewportShowExpand: false,
    viewportBackgroundColor: background,
  })
  viewer.plugin.canvas3d?.setProps({
    // the camera is aimed from Python; Mol*'s own reset on load would undo it
    camera: { helper: { axes: { name: "off", params: {} } }, manualReset: true,
      mode: "orthographic" },
    transparentBackground: true,
    checkeredTransparentBackground: false,
  })
}

/** color: a number (uniform), or {theme, params} for any Mol* colour theme;
 *  size: optional {theme, params} for the size theme. */
async function addRepr(plugin, structure, expression, key, color, type, typeParams, size) {
  const component = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure, expression, key, { label: key })
  if (!component) return null
  const theme = typeof color === "number"
    ? { color: "uniform", colorParams: { value: color } }
    : { color: color.theme, colorParams: color.params }
  await plugin.builders.structure.representation.addRepresentation(component, {
    type, typeParams, ...theme,
    ...(size ? { size: size.theme, sizeParams: size.params } : {}),
  })
  return component
}

const addBallAndStick = (plugin, structure, expression, key, color, size, alpha = 1) =>
  addRepr(plugin, structure, expression, key, color, "ball-and-stick", { sizeFactor: size, alpha })

/** Sticks with no balls: every atom the radius of its bonds. Mol*'s default
 *  size theme gives each atom its van der Waals radius, which leaves knobs on
 *  O, N and S; a uniform size of STICK_R (about a carbon) keeps the bonds as
 *  thick as before and the joints flush. Aromatic rings and double bonds are
 *  drawn as single sticks too. */
const STICK_R = 1.7
const STICK_SIZE = 0.2                  // one thickness for catalytic, secondary and heme; colour tells them apart
const addSticks = (plugin, structure, expression, key, color, size, alpha = 1) =>
  addRepr(plugin, structure, expression, key, color, "ball-and-stick",
    { sizeFactor: size, sizeAspectRatio: 1, alpha, aromaticBonds: false, multipleBonds: "off" },
    { theme: "uniform", params: { value: STICK_R } })

/**
 * spec: {pdb (text), keep [{chain, beg, end}], sites [{chain, num, kind, ok,
 * label, anchor [x,y,z]}], cofactors [{chain, num}], ligand bool, color,
 * face [x,y,z], up [x,y,z], centre [x,y,z]}
 * view: {zoom bool, size css px, zoomRadius angstrom}
 */
window.renderScene = async (spec, view) => {
  const { Schema } = molstar.lib.structure.StructureElement
  const { Structure } = molstar.lib.structure
  if (!viewer) await boot(view.background || "#ffffff")
  const plugin = viewer.plugin
  await plugin.clear()
  // optional, for colour comparisons (the mechanism-groups demo): spec.palette
  // overrides COLORS; spec.secondary_size scales the secondary sticks;
  // spec.substrate_style "sticks" draws the plastic as thin sticks;
  // spec.outline draws Mol*'s dark outline around everything
  const C = { ...COLORS, ...Object.fromEntries(Object.entries(spec.palette || {}).map(([k, v]) => [k, hex(v)])) }
  if (spec.outline) plugin.canvas3d?.setProps({ postprocessing: { outline:
    { name: "on", params: { scale: 1, threshold: 0.33, color: 0x2b2f36,
      includeTransparent: !spec.see_through } } } })   // see-through parts stay unoutlined

  const data = await plugin.builders.data.rawData({ data: spec.pdb, label: spec.label })
  const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb")
  const model = await plugin.builders.structure.createModel(trajectory)
  const structure = await plugin.builders.structure.createStructure(model)
  const sel = (items) => Schema.toExpression({ items })

  // the folded core only: dangling predicted termini are cut on the Python
  // side (the step 12 rule). One cartoon per protein chain, so a partner
  // chain can be told apart
  // the residues in front of the mechanism are cut away, in both views
  // (zoom_keep is the drawn chain minus those)
  const colorOf = Object.fromEntries(spec.keep.map((k) => [k.chain, k.color]))
  const pieces = spec.zoom_keep || spec.keep
  for (const chain of Object.keys(colorOf)) {
    const items = pieces.filter((k) => k.chain === chain)
      .map((k) => ({ auth_asym_id: k.chain, beg_auth_seq_id: k.beg, end_auth_seq_id: k.end }))
    if (!items.length) continue
    await addRepr(plugin, structure, sel(items), `chain ${chain}`,
      hex(spec.cartoon_override || colorOf[chain]),
      "cartoon", { alpha: view.zoom ? 0.6 : 1 })
  }
  // optional: backbone in front of the mechanism, see-through instead of cut away
  for (const chain of Object.keys(colorOf)) {
    const items = (spec.see_through || []).filter((k) => k.chain === chain)
      .map((k) => ({ auth_asym_id: k.chain, beg_auth_seq_id: k.beg, end_auth_seq_id: k.end }))
    if (!items.length) continue
    await addRepr(plugin, structure, sel(items), `see-through ${chain}`,
      hex(spec.cartoon_override || colorOf[chain]), "cartoon", { alpha: spec.see_through_alpha })
  }

  const big = view.zoom ? 1 : 1.6        // sticks read at thumbnail size too
  const group = (kind, ok) => spec.sites.filter((s) => s.kind === kind && s.ok === ok)
  const pick = (s) => sel(s.map((r) => ({ auth_asym_id: r.chain, auth_seq_id: r.num })))
  for (const ok of [true, false]) {
    const sec = group("secondary", ok), cat = group("primary", ok)
    if (sec.length) await addSticks(plugin, structure, pick(sec), `secondary${ok ? "" : "?"}`,
      C.secondary, STICK_SIZE * big * (spec.secondary_size || 1), ok ? 1 : 0.4)
    if (cat.length) await addSticks(plugin, structure, pick(cat), `catalytic${ok ? "" : "?"}`,
      C.catalytic, STICK_SIZE * big * (spec.catalytic_size || 1), ok ? 1 : 0.4)
  }
  const metals = spec.cofactors.filter((c) => c.metal), organic = spec.cofactors.filter((c) => !c.metal)
  if (organic.length) await addSticks(plugin, structure, pick(organic), "heme",
    C.heme, STICK_SIZE * big)
  if (metals.length) await addRepr(plugin, structure, pick(metals), "metal ions",
    C.metal, "spacefill", { sizeFactor: 0.75 })
  if (spec.ligand) {
    const L = sel([{ auth_asym_id: "L" }])
    // substrate_carbon: element colours (O red, N blue) with carbons in that
    // colour; substrate_color: one colour. Both optional, for comparisons
    const sub = spec.substrate_carbon
      ? { theme: "element-symbol", params: { carbonColor: { name: "uniform",
          params: { value: hex(spec.substrate_carbon) } } } }
      : hex(spec.substrate_color || C.substrate)
    const front = spec.substrate_front || []
    if (spec.substrate_style === "sticks" && front.length) {
      // the whole plastic see-through, then its atoms behind the mechanism solid
      // over it, so the bonds across the split are still drawn
      await addSticks(plugin, structure, L, "plastic (see-through)", sub, STICK_SIZE * big * 0.7,
        spec.substrate_alpha ?? spec.see_through_alpha)
      const back = sel((spec.substrate_back || []).map((n) => ({ auth_asym_id: "L", auth_atom_id: n })))
      if ((spec.substrate_back || []).length)
        await addSticks(plugin, structure, back, "plastic", sub, STICK_SIZE * big * 0.7)
    } else if (spec.substrate_style === "sticks")
      await addSticks(plugin, structure, L, "plastic", sub, STICK_SIZE * big * 0.7)
    else
      await addRepr(plugin, structure, L, "plastic", sub, "spacefill", { sizeFactor: 0.7 })
  }

  // frame: the core's box across the screen, or the site and the plastic,
  // looking along -face so the active site faces the reader
  const camera = plugin.canvas3d.camera
  const dir = spec.face.map((v) => -v)
  const r = view.zoom ? view.zoomRadius : spec.frame_radius * (spec.frame_scale || 1)
  const at = view.zoom ? (spec.zoom_centre || spec.centre) : spec.frame_centre
  const focus = camera.getInvariantFocus(at, r, spec.up, dir)
  // only these fields: the rest of the snapshot belongs to the last structure
  // radius sets the clipping planes, not the framing (that is the camera
  // distance), so keep it wide enough that nothing drawn is sliced
  const state = {
    target: focus.target, position: focus.position, up: focus.up,
    radius: view.zoom ? spec.scene_radius : r,
    radiusMax: spec.scene_radius,
  }
  await wait(300)
  camera.setState(state, 0)
  plugin.canvas3d.requestDraw(true)
  // occlusion and anti-aliasing refine over a few frames
  await wait(1200)

  // project points for the overlay: canvas pixels, origin bottom left
  const canvas = plugin.canvas3d.webgl.gl.canvas
  const box = canvas.getBoundingClientRect()
  const toCss = (p) => {
    const out = new Float32Array(4)
    camera.project(out, new Float32Array(p))
    return [out[0] / canvas.width * box.width, box.height - out[1] / canvas.height * box.height]
  }
  const zoomAt = spec.zoom_centre || spec.centre
  const centre = toCss(zoomAt)
  const ring = (() => {
    // the zoom window as seen in the full view, from a point zoomRadius off
    // the centre across the screen
    const right = cross(spec.face, spec.up)
    const edge = toCss(zoomAt.map((c, i) => c + right[i] * view.zoomRadius))
    return { x: centre[0], y: centre[1], r: Math.hypot(edge[0] - centre[0], edge[1] - centre[1]) }
  })()
  const labels = []
  if (view.labels) {
    const seen = new Set()
    const inView = ([x, y]) => x > 10 && y > 10 && x < view.size - 10 && y < view.size - 10
    for (const s of spec.sites) {
      // a residue listed in two groups is labelled once; one outside the
      // window is not labelled at all
      if (seen.has(s.label) || !inView(toCss(s.anchor))) continue
      seen.add(s.label)
      const [x, y] = toCss(s.anchor)
      let dx = x - centre[0], dy = y - centre[1]
      const n = Math.hypot(dx, dy) || 1
      dx /= n; dy /= n
      labels.push({ x: x + dx * 16, y: y + dy * 16, ax: x, ay: y, text: s.label + (s.ok ? "" : "?"), kind: s.kind })
    }
    // the plastic and the cofactors, named where they sit (inside the window)
    const inside = ([x, y]) => x > 20 && y > 20 && x < view.size - 20 && y < view.size - 20
    if (spec.ligand_anchor) {
      const p = toCss(spec.ligand_anchor)
      if (inside(p)) labels.push({ x: p[0], y: p[1] + 30, text: spec.ligand_short || spec.ligand_name, kind: "substrate" })
    }
    // a cofactor much bigger than a residue (heme) is named at its rim, on the
    // side away from the labelled residues, so its name does not sit on theirs
    const prim = spec.sites.map((q) => toCss(q.anchor)).filter(inView)
    const pc = prim.length ? [prim.reduce((a, q) => a + q[0], 0) / prim.length,
      prim.reduce((a, q) => a + q[1], 0) / prim.length] : null
    const named = new Set()
    for (const c of spec.cofactors) {
      const p = toCss(c.anchor)
      const key = `${c.name}|${Math.round(p[0] / 40)}|${Math.round(p[1] / 40)}`
      if (!inside(p) || named.has(key)) continue
      named.add(key)
      let ox = 26, oy = -22
      if (!c.metal && pc) {
        let ex = p[0] - pc[0], ey = p[1] - pc[1]
        const en = Math.hypot(ex, ey)
        if (en > 1) { ox = ex / en * 150; oy = ey / en * 150 } else { ox = 0; oy = -150 }
      }
      labels.push({ x: p[0] + ox, y: p[1] + oy, text: c.metal ? c.name.replace(" ion", "") : c.name, kind: "cofactor" })
    }
    spreadLabels(labels, view.size)
    // hideLabels: positions only (views.json), for labels drawn elsewhere (step 18 insets)
    drawLabels(view.hideLabels ? [] : labels)
  } else {
    drawLabels([])
  }
  return { ring, labels, size: view.size }
}

function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/** Push overlapping labels apart, a few relaxation rounds, kept on canvas. */
function spreadLabels(labels, size) {
  const h = 26
  // width from the text, so a long name is kept inside the frame too
  for (const l of labels) l.w = Math.max(64, l.text.length * 14)
  for (let it = 0; it < 60; it++) {
    for (let i = 0; i < labels.length; i++) {
      for (let j = i + 1; j < labels.length; j++) {
        const a = labels[i], b = labels[j]
        const ox = (a.w + b.w) / 2 - Math.abs(a.x - b.x), oy = h - Math.abs(a.y - b.y)
        if (ox > 0 && oy > 0) {
          if (oy < ox) {
            const s = (a.y < b.y ? -1 : 1) * oy / 2
            a.y += s; b.y -= s
          } else {
            const s = (a.x < b.x ? -1 : 1) * ox / 2
            a.x += s; b.x -= s
          }
        }
      }
    }
    for (const l of labels) {
      l.x = Math.min(size - l.w / 2 - 6, Math.max(l.w / 2 + 6, l.x))
      l.y = Math.min(size - h / 2, Math.max(h / 2, l.y))
    }
  }
}

function drawLabels(labels) {
  const layer = document.getElementById("labels")
  layer.innerHTML = ""
  for (const l of labels) {
    const el = document.createElement("div")
    el.className = `label ${l.kind}`
    el.textContent = l.text
    el.style.left = `${l.x}px`
    el.style.top = `${l.y}px`
    layer.appendChild(el)
  }
}
