# Cómo publicarlo / How to publish

Repo name assumed throughout: **`plasmonics-multilayer`** under **`KANKALU`**.
If you rename it, update the three links in `docs/index.html` and the two in `README.md`.

## 1. Create the repository

On GitHub: **New repository** → name `plasmonics-multilayer` → **Public** →
do *not* add a README, .gitignore or licence (they are already here).

## 2. Push

```bash
cd plasmonics-multilayer
git init
git add .
git commit -m "Transfer-matrix reflectance of metal-oxide multilayers"
git branch -M main
git remote add origin https://github.com/KANKALU/plasmonics-multilayer.git
git push -u origin main
```

## 3. Turn on GitHub Pages

Repository → **Settings** → **Pages** → Source: *Deploy from a branch* →
Branch: `main`, folder: **`/docs`** → **Save**.

A minute later the site is at:

    https://kankalu.github.io/plasmonics-multilayer/

## 4. Add the link to the repo header

Repository home → the gear next to **About** → paste the Pages URL into *Website*,
and add topics: `physics`, `optics`, `plasmonics`, `thin-films`, `transfer-matrix`,
`numpy`, `scientific-computing`.

## Previewing the site locally

`fetch()` will not read `data/curves.json` from a `file://` URL, so open it over HTTP:

```bash
cd docs && python -m http.server 8000
# then visit http://localhost:8000
```

## If you change the physics code

`docs/data/curves.json` is generated, not hand-written. After editing anything under
`src/plasmonics/`, regenerate it so the site keeps matching the library:

```bash
python scripts/export_web_data.py
```

## Not included on purpose

- **The raw goniometer files** (~400 MB of `*_Norm_Air.txt` / `.csv`). `.gitignore` blocks
  them by pattern, so they cannot be committed by accident. If you ever want them public,
  attach them to a GitHub **Release** rather than tracking them in the repo.
- **`PLASMONICS_Maier_Springer.pdf`.** That book is copyrighted; it is cited in the README
  and in both reports, and must not be uploaded.
