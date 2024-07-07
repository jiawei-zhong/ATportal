<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Single cell.png"/></div>
    <div class="column-2-mod">
        <div class="mod_intro_text">
            <p>This module contains single-cell and single-nucleus RNA sequencing results from visceral and subcutaneous WAT. In total, data from XX cross-sectional studies, which together compriseing ZZZZ individuals (YY men, XX women), are available.
            </p>
          </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0"  src="https://atportal.medh.ki.se/shiny/Module_Singlecell/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModuleSingleCell'
  }

  </script>
  
  <style scoped>

  </style>
  