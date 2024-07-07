<template>
    <div>
      <HeaderComponent />
      <div class="app-container">
        <main>

<div class="banner_mod">
  <div class="grid-container_mod">
    <div class="column-1-mod"><img src="@/assets/module_images_Perturbations.png"/></div>
    <div class="column-2-mod">
        <div class="mod_intro_text">
            <p>This module contains transcriptomic data detailing the effects of different perturbations on gene expression in human adipocytes.<br>
    The following user-defined queries can be interrogated in this module:<br>
    <b>External perturbations: </b>explore transcript levels in adipocytes following external perturbations such as hypoxia and TNF-a incubations.<br>
    <b>Genetic perturbations: </b>study transcriptional changes in adipocytes following depletion of individual genes.
            </p>
          </div>
    </div>
  </div>
</div>
<hr class="separation_mod">


<div class="iframe_pos">
  <iframe height="1000" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Perturbation_v03/"></iframe>
</div>

</main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  export default {
    name: 'ModulePerturbation'
  }

  </script>
  
  <style scoped>

  </style>
  