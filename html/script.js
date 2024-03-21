document.addEventListener("DOMContentLoaded", function () {
    const stickyColumn = document.getElementById("sticky-column");
    const footer = document.getElementById("footer");

    window.addEventListener("scroll", function () {
        const scrollTop = window.scrollY;
        const maxTranslation = footer.offsetTop - stickyColumn.offsetHeight;
        const translateY = Math.min(scrollTop, maxTranslation);

        stickyColumn.style.transition = "transform 0.005s ease";
        stickyColumn.style.transform = `translateY(${translateY}px)`;
    });
});