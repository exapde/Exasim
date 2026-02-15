# Exasim Showcase Website

This folder contains a GitHub Pages-ready static site for showcasing Exasim simulation images and videos.

## Local preview

From the repository root:

```bash
cd website
python3 -m http.server 8000
```

Then open `http://localhost:8000`.

## Add images

1. Copy images into `website/assets/images/`.
2. Add entries in `website/content.js` under `images`.

## Add videos

1. Copy `.mp4` or `.webm` files into `website/assets/videos/`.
2. Add entries in `website/content.js` under `videos` with `src: "assets/videos/<file>.mp4"`.

## Deploy on GitHub Pages

1. Push changes to `master`.
2. In GitHub repo settings, open **Pages**.
3. Set **Build and deployment** source to **GitHub Actions**.
4. The workflow in `.github/workflows/exasim-showcase-pages.yml` deploys the `website/` folder.

Your site URL will be:

- `https://exapde.github.io/Exasim/`
