<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Spatial.png"/></div>
    <div class="column-2-mod">
        <div class="mod_intro_text">
            <p>This module contains spatial transcriptomic data detailing the cytoarchitecture of subcutaneous WAT. In total, data from one cross-sectional study comprising 10 individuals (YY men, XX women), are available.
            </p>
          </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Spatial/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModuleSpatial'
  }

  </script>
  
  <style scoped>

  </style>
  