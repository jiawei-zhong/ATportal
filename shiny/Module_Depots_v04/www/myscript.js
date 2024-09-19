const sidebar = document.querySelector('.sidebar');
const mainContent = document.querySelector('.main-content');
document.querySelector('.trigger').onclick = function () {
  sidebar.classList.toggle('sidebar_small');
  mainContent.classList.toggle('main-content_large')
}