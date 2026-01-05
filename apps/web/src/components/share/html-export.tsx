import type { Atom } from '../../lib/types'
import { atomsToPdb } from '../../lib/pdb'

function buildHtml(pdbText: string): string {
  const safePdb = pdbText.replace(/`/g, '\\`')
  return `<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Chem Model Share</title>
  <script src="https://unpkg.com/ngl@2/dist/ngl.js"></script>
  <style>
    html, body { margin: 0; height: 100%; background: #0b1120; }
    #viewport { width: 100vw; height: 100vh; }
    .badge { position: absolute; top: 16px; left: 16px; padding: 6px 12px; background: rgba(15,23,42,0.7); color: #e2e8f0; border-radius: 999px; font-family: system-ui, sans-serif; font-size: 12px; }
  </style>
</head>
<body>
  <div class="badge">Chem Model Share (NGL)</div>
  <div id="viewport"></div>
  <script>
    const stage = new NGL.Stage('viewport', { backgroundColor: '#0b1120' });
    const pdbData = ` + "`" + `${safePdb}` + "`" + `;
    const blob = new Blob([pdbData], { type: 'text/plain' });
    stage.loadFile(blob, { ext: 'pdb', name: 'Model' }).then((comp) => {
      comp.addRepresentation('ball+stick', { opacity: 1.0 });
      stage.autoView();
    });
  </script>
</body>
</html>`
}

export function downloadShareHtml(atoms: Atom[], filename = 'chem-model-share.html') {
  const pdb = atomsToPdb(atoms)
  const html = buildHtml(pdb)
  const blob = new Blob([html], { type: 'text/html' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  URL.revokeObjectURL(url)
}
