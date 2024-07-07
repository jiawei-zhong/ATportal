<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Clinical.png"/></div>
    <div class="column-2-mod">
      <div class="mod_intro_text">
        <p>This module contains clinical/adipocyte parameters and WAT transcriptome/proteome results from subcutaneous and omental WAT from XX cross-sectional and YY longitudinal studies together comprising ZZZZ individuals (YY men, XX women). 
        To learn more about the different cohorts click <a href="vignettes.html"
            target="_blank">here</a>. <br>
          The following user-defined queries can be interrogated in this module: <br>
          <b>Clinical associations:</b> identify associations to clinical and adipocyte parameters.<br>
          <b>Sex differences:</b> explore differences between women and men.<br>
          <b>Obesity:</b> differences between individuals with or without obesity.<br>
          <b>Weight loss:</b> define effects of short- and/or long-term weight loss.<br>
        </p>
      </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Clinical_v03/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModuleClinical'
  }

  </script>
  
  <style scoped>


  </style>
  