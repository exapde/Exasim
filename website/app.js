(function () {
  const content = window.exasimContent || { images: [], videos: [] };
  const imageGrid = document.getElementById("image-grid");
  const videoGrid = document.getElementById("video-grid");
  const filtersWrap = document.getElementById("filters");

  const categories = ["All", ...new Set(content.images.map((item) => item.category))];
  let activeCategory = "All";

  function cardTemplate(item) {
    return `
      <article class="card reveal" style="--delay: 0.05s;">
        <img src="${item.src}" alt="${item.title}">
        <div class="card-body">
          <h3>${item.title}</h3>
          <p>${item.description}</p>
          <span class="tag">${item.category}</span>
        </div>
      </article>
    `;
  }

  function videoTemplate(item) {
    const media = item.src
      ? `<video controls preload="metadata" poster="${item.poster || ""}"><source src="${item.src}"></video>`
      : `<img src="${item.poster || "assets/images/logo.png"}" alt="${item.title}">`;

    return `
      <article class="card reveal" style="--delay: 0.05s;">
        ${media}
        <div class="card-body">
          <h3>${item.title}</h3>
          <p>${item.description}</p>
          <span class="tag">${item.category}</span>
        </div>
      </article>
    `;
  }

  function renderFilters() {
    filtersWrap.innerHTML = categories
      .map(
        (category) =>
          `<button class="filter ${category === activeCategory ? "is-active" : ""}" data-category="${category}">${category}</button>`
      )
      .join("");

    filtersWrap.querySelectorAll(".filter").forEach((button) => {
      button.addEventListener("click", () => {
        activeCategory = button.dataset.category;
        renderGallery();
        renderFilters();
      });
    });
  }

  function renderGallery() {
    const filtered =
      activeCategory === "All"
        ? content.images
        : content.images.filter((item) => item.category === activeCategory);

    imageGrid.innerHTML = filtered.map(cardTemplate).join("");
    observeReveal();
  }

  function renderVideos() {
    videoGrid.innerHTML = content.videos.map(videoTemplate).join("");
    observeReveal();
  }

  function observeReveal() {
    const targets = document.querySelectorAll(".reveal");
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add("in-view");
            observer.unobserve(entry.target);
          }
        });
      },
      { threshold: 0.12 }
    );

    targets.forEach((target) => observer.observe(target));
  }

  renderFilters();
  renderGallery();
  renderVideos();
  observeReveal();
})();
