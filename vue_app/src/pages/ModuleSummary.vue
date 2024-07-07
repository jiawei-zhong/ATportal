<template>
    <div>
      <HeaderComponent />
      <BannerComponent />
      <div class="app-container">
        <main>

            <div class="grid-container">
      <div class="column-1" id="sticky-column">
        <!-- <div class="center-container2"> -->
         <!-- <div class="search-container2"> -->
            <!-- <input type="text" id="userInput" class="search-bar2" placeholder="Search gene"> -->
            <!-- <button onclick="sendToShiny()" class="search-button2">Search</button> -->
          <!-- </div> -->
        <!-- </div> -->
        <!--<hr class="separation_sum_side">-->
        <div class="text-container">
          <p2 id="gene-description"></p2>
        </div>
        <div class="download-container">
          <button id="download-data" class="download-button2" style="color: black; max-width: 200px;">Download</button>
          <!-- <button id="download-data" class="btn btn-primary" style="color: black; max-width: 200px;" onclick="downloadPdf()">Download</button> -->
        </div>
      </div>

      <div class="column-2">
        <div class="iframe-group">

          <div class="iframe-title">Clinical correlation
            <a href="clinical.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe1" src="https://atportal.medh.ki.se/shiny/Summary_gene_Clinical_Correlation/"
            height="380" class="collapsible-iframe"></iframe>

          <div class="iframe-title">Sex difference
            <a href="clinical.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe2" src="https://atportal.medh.ki.se/shiny/Summary_gene_Clinical_Sex/" height="400"
            class="collapsible-iframe"></iframe>

          <div class="iframe-title">Obese vs. non-obese
            <a href="clinical.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe3" src="https://atportal.medh.ki.se/shiny/Summary_gene_Clinical_Obese/" height="400"
            class="collapsible-iframe"></iframe>

          <div class="iframe-title">Weight loss
            <a href="clinical.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe4" src="https://atportal.medh.ki.se/shiny/Summary_gene_Weightloss/" height="400"
            class="collapsible-iframe"></iframe>
            
          <div class="iframe-title">Depots
            <a href="depots.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe5" src="https://atportal.medh.ki.se/shiny/Summary_gene_Depots/" height="250"
            class="collapsible-iframe"></iframe>

          <div class="iframe-title">Adipogenesis
            <a href="celltype.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe6" src="https://atportal.medh.ki.se/shiny/Summary_gene_Adipogenesis/" height="250"
            class="collapsible-iframe"></iframe>

        <div class="iframe-title">Tissue specificity
             <a href="spatial.html" target="_blank">
             <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
              </a>
          </div>
          <iframe id="shinyIframe7" src="https://atportal.medh.ki.se/shiny/Summary_gene_Tissue/" height="250"
            class="collapsible-iframe"></iframe>

          <div class="iframe-title">Single cell
            <a href="sc.html" target="_blank">
              <ion-icon name="caret-forward-circle" class="arrow-icon"></ion-icon>
            </a>
          </div>
          <iframe id="shinyIframe8" src="https://atportal.medh.ki.se/shiny/Summary_gene_Singlecell/" height="420"
            class="collapsible-iframe"></iframe>

        </div>

      </div>
    </div>

        </main>
      </div>
      <FooterComponent />
    </div>
  </template>
  
  <script>
  import BannerComponent from '@/components/BannerComponent.vue';

  export default {
    name: 'ModuleSummary',
     components: {
        BannerComponent
     }
  }

  </script>
  
  <style scoped>


.arrow-icon {
    color: #1a4659;
    vertical-align: middle;
}   

.grid-container {
    padding-top: 2%;
    margin: auto;
    display: flex;
    /* overflow: hidden; */
    /* position: relative; */
    /* justify-content: space-between; */
    width: 81%;
    height: 100%;
    grid-template-columns: min-content 1fr;
    grid-gap: 1px;
}


.column-1 {
    flex: 1;
    flex-direction: column;
    position: sticky;
    top: 0; /* Make the sidebar stick to the top of the viewport */
    display: grid;
    /* max-height: 40rem; */
    height:100%;
    grid-template-rows: auto 1fr;
    background-color: #edf2f2;
    padding: 1rem;
    border:0.15rem solid rgba(26,70,89,1.00);
    border-radius: 2rem;
    /* overflow-y: auto;  */
    overflow: auto;
    transition: transform 0.005s ease;
}

.center-container2 {
    display: flex;
    justify-content: center;
    align-items: center;
    padding-bottom: 0.9rem;
    padding-top: 0.5rem;
}

.search-container2 {
    display: flex;
    align-items: center;
    text-align: center;
    padding: 0rem;
    width: 90%;
    margin: auto;
}

.search-bar2 {
    background: rgba(26, 70, 89, 0.8);
    color: #edf2f2;
    width: 120%;
    border: 0.03rem solid rgba(26,70,89,1.00);
    border-radius: 2rem;
    padding: 0.6rem;
    outline: none;
}

.search-button2 {
    background-color: #e2c744;
    width: 25%;
    color: #edf2f2;
    border: none;
    border-radius: 2rem;
    cursor: pointer;
    padding: 0.6rem;
    margin-left: 0.5rem;
}


.search-button2 &::placeholder {
    color: #1a4659;
    text-align: justify;
}

.separation_sum_side{
    width: 90%;
    max-height: 0.5px;
    background-color: #1a4659;
    border: 0.01rem solid rgba(26,70,89,1.00);
}

.text-container {
    width: 90%;
    height: auto;
    padding-left: 0.6rem;
    padding-right: 0.6rem;
    text-align: justify-all;
    color: #1a4659;
    font-size: 1rem;
}

p2 {font-size: 0.9rem}

.download-container2 {
    display: flex;
    align-items: center;
    text-align: center;
    padding: 0rem;
    width: 90%;
    margin: auto;
}

.download-button2 {
    background-color: #e2c744;
    width: 80px;
    color: #edf2f2 !important; /* 强制应用颜色 */
    border: none;
    border-radius: 2rem;
    cursor: pointer;
    padding: 0.6rem;
    margin-left: 0.5rem;
}


.download-button2 &::placeholder {
    color: #c7d0d3;
    text-align: justify;
}



.column-2 {
    flex: 3;
    overflow-x: auto;  /* 当内容宽度超出时显示水平滚动条 */
    flex-direction: column;
    padding-left: 2rem;
    padding-top: 0rem;
    height: 100%;
}

/*
#sticky-hr {
    margin-top: 0px;
  }
*/

  </style>
  