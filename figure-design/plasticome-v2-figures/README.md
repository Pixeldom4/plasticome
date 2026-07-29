# Plasticome V2 Rebuild — pipeline figures

Flowchart of the V2 rebuild (Step 1 → Step 2.6), rendered in four diagram languages.
Same pipeline in each; pick the one whose look/tooling you prefer.

| Language | Source | Raster | Vector | Notes |
|----------|--------|--------|--------|-------|
| Mermaid  | `pipeline.mmd`  | `mermaid.png`  | `mermaid.svg`  | Renders natively in Obsidian (no plugin). Most detail: seed/no-seed split + 2.6a/b/c. |
| D2       | `pipeline.d2`   | `d2.png`       | `d2.svg`       | Cleanest export for slides/papers. Needs `d2` CLI or D2 plugin. |
| Graphviz | `pipeline.dot`  | `graphviz.png` | `graphviz.svg` | Most compact/dense. Needs `dot` (graphviz) or plugin. |
| PlantUML | `pipeline.puml` | `plantuml.png` | `plantuml.svg` | Linear activity style; does not branch the two input arms. Needs Java + plantuml.jar. |

## Re-render

```sh
d2 --theme 0 pipeline.d2 d2.svg
dot -Tsvg pipeline.dot -o graphviz.svg
java -jar plantuml.jar -tsvg pipeline.puml
npx @mermaid-js/mermaid-cli -i pipeline.mmd -o mermaid.svg -b white
```

Recommendation: **Mermaid** for the vault (zero-setup, correct Y-merge), **D2** for a polished exported figure.
