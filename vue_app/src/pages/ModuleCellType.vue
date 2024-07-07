<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Cell type.png"/></div>
    <div class="column-2-mod">
        <div class="mod_intro_text">
            <p>This module contains transcriptomic data detailing expression in different cell types within WAT and other tissues as well as during adipocyte differentiation. 
			To learn more about the different cohorts click <a href="vignettes.html" target="_blank">here</a>. <br>
              The following user-defined queries can be interrogated in this module: <br>
              <b>Tissue comparisons:</b> explore expression patterns across multiple human tissues.<br>
              <b>Adipogenesis</b>investigate expression levels during adipocyte formation. <br>
            </p>
          </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Characterization_v06/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModuleCellType'
  }

  </script>
  
  <style scoped>

  </style>
  