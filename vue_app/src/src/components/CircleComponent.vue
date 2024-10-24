<template>
 <div id="main-container">

<div class="circles">

  <div id="modules1">
    <div class="grid-container_image">
      <div class="grid-child">
        <router-link class="tooltip" to="/summary"><img src="@/assets/images/module_images_Summary.png"
            alt="Summary"><span>Obtain summarized information from all aspects of the portal.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/clinical">
          <img src="@/assets/images/module_images_Clinical.png" alt="Clinical"><span>Identify associations to clinical and
            fat cell parameters.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/depots">
            <img src="@/assets/images/module_images_Depots.png"
            alt="Depots"><span>Compare subcutaneous and visceral WAT depots.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/celltype">
            <img src="@/assets/images/module_images_Cell type.png"
            alt="Cell_type"><span>Define enrichment in multiple cell populations.</span></router-link>
      </div>
    </div>
  </div>


  <div id="modules2">
    <div class="grid-container_image">
      <div class="grid-child">
        <router-link class="tooltip" to="/singlecell">
            <img src="@/assets/images/module_images_Single cell.png"
            alt="Single_cell"><span>Explore single-cell and single-nucleus RNAseq. results.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/spatial">
            <img src="@/assets/images/module_images_Spatial.png"
            alt="Spatial"><span>Investigate spatial transcriptomics data.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/perturbation">
            <img src="@/assets/images/module_images_Perturbations.png"
            alt="Perturbation"><span>Determine the influence of different perturbations.</span></router-link>
      </div>
      <div class="grid-child">
        <router-link class="tooltip" to="/genefinder"><img src="@/assets/images/module_images_Gene_finder.png"
            alt="Gene_finder"><span>Discover genes/proteins fulfilling user-defined criteria.</span></router-link>
      </div>
    </div>
  </div>

</div>
</div>
  </template>
  
  <script>
  export default {
    name: 'CircleComponent'
  }
  </script>
  
  <style scoped>
  #main-container {
    display: table;
    width: 100%;
    height: 100%;
    margin-top: 1rem;
    text-align: center;
}
.circles{
    width: 80%;
    max-width: 1400px;
    min-width: 500px;
    margin: auto;
    text-align: center;
    align-self: center;
}

.grid-container_image {
    width: 100%;
    display: grid;
    grid-template-columns: 1fr 1fr 1fr 1fr;
    grid-gap: 0px;
    margin: auto;
    padding-left: 0px;
    padding-right: 0px;
    position: relative;
    text-align: center;
}

.grid-child a img {
    width: 80%;
    max-width: 480px;
}

.grid-child a text {
    width: 70%;
    max-width: 480px;
}

a.tooltip {
    border-bottom: px;
    text-decoration: none;
}
a.tooltip:hover {
    cursor: help;
    position:inherit;
}
a.tooltip span {
    display: none;
}
a.tooltip:hover span {
    border: #1a4659 0.125rem;
    border-radius: 1.125rem;
    padding: 1% 1% 1% 1%;
    display: block;
    transition: transform 250ms ease-in-out;
    background: #edf2f2;
    margin: 0rem;
    width: 25%;
    position: absolute;
    top: 50%;
    text-decoration: none;
}



  </style>