<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Depots.png"/></div>
    <div class="column-2-mod">
        <div class="mod_intro_text">
            <p>This module allows for comparisons of gene and protein expression in visceral and subcutaneous WAT depots. In total, data from XX cross-sectional studies, which together compriseing ZZZZ individuals (YY men, XX women), are available.
			To learn more about the different cohorts click <a href="vignettes.html" target="_blank">here</a>. <br></p>
          </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Depots_v02/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModuleDepots'
  }

  </script>
  
  <style scoped>

  </style>
  