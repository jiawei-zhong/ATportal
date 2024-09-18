<template>
    <div class="header_col">
      <div class="banner">
        <div><router-link to="/"><img src="@/assets/images/logo_banner_full_2.png" alt="" class="logo_large" /></router-link></div>
      </div>
      <hr class="separation">
      <div class="intro_text">Open-access to clinical, transcriptional and protein data down to the single-cell level</div>
      <form id="myForm" class="search_bar" @submit.prevent="modifyAction">
        <div class="center-container">
          <div class="search-container">
            <div class="search-bar-container">
              <div class="search-icon">
                <!-- SVG icon here -->
              </div>
              <input type="text" id="gene" name="gene" v-model="searchGene" class="search-bar" placeholder="Search gene, eg. LEP">
            </div>
            <button type="submit" class="search-button">Search</button>
          </div>
        </div>
      </form>
    </div>
  </template>
  
  <script>
  import $ from 'jquery';
  import 'jquery-ui/ui/widgets/autocomplete';
  
  export default {
    name: 'BannerComponent',
    data() {
      return {
        searchGene: ''
      };
    },
    mounted() {
      this.$nextTick(() => {
        this.initAutocomplete();
      });
    },
    methods: {
      initAutocomplete() {
        $("#gene").autocomplete({
          source: (request, response) => {
            $.ajax({
              url: "/autocomplete",
              dataType: "json",
              data: {
                term: request.term
              },
              success: data => response(data)
            });
          },
          minLength: 2,
          select: (event, ui) => {
            this.searchGene = ui.item.value;
            this.modifyAction();
            return false;
          }
        }).autocomplete("instance")._resizeMenu = function() {
          this.menu.element.outerWidth(this.element.outerWidth());
        };
      },
      modifyAction() {
        const targetUrl = `/summary?gene=${encodeURIComponent(this.searchGene)}`;
        this.$router.push(targetUrl);
      }
    }
  }
  </script>
  
  
  
  <style scoped>
  .banner * {
    display: inline-block;
    vertical-align: middle;
    padding-left: auto;
    padding-right: auto;
    margin-left: px;
    margin-right: auto;
    margin-top: 0.5rem;
}

.header_col {
    background-color: #1a4659;
    width: auto;
    min-width: 500px;
}


.intro_text{
    margin: auto;
    text-align: center;
    font-size: 124%;
    padding-bottom: 0.5rem;
    color:#edf2f2;
}

.center-container {
    display: flex;
    justify-content: center;
    align-items: center;
    padding-bottom: 0.1rem;
}

.search-container {
    display: flex;
    align-items: center;
    /* border: px solid #ccc; */
    /* border-radius: 5px; */
    padding: 1%;
    width: 43.7rem;
}

.search-bar-container {
    display: flex;
    align-items: center;
    flex: 1;
    border: 0.01rem solid #edf2f2;
    border-radius: 3rem;
    background: rgba(237, 242,242, 0.1);
}

.search-bar {
    border: none;
    width: 100%;
    padding: 1%;
    outline: none;
    background: transparent;
    color: #edf2f2;
}

.search-icon {
    display: flex;
    align-items: center;
    padding: 1%;
}

.search-icon svg {
    width: 2rem;
    height: 1.2rem;
    fill: #edf2f2;
}

.search-button {
    background-color: #e2c744;
    font: "Red Hat Display";
    font-weight: 300;
    color: #edf2f2;
    border: none;
    border-radius: 3rem;
    cursor: pointer;
    width: 12%;
    padding: 1.2%;
    margin-left: 1.2%;
}

.search-button:hover {
    background-color: rgba(237, 242, 242, 0.8);
}

::placeholder {
    color: white;
    text-align: justify;
}


.banner    {
    padding-top: rem;
    text-align: center;
    background-color: transparent;
}

.banner_text {
    width: 26.25rem;
    margin-top: 1.875rem;
    padding-left: 0.625rem;
    text-align: left;
    font-size: 2rem;
}

.logo_large {
    width: 32rem;
    right: auto;
    left: auto;
}

.separation {
    width: 43.7rem;
    margin-top: 1.25rem;
    margin-bottom: 1.25rem;
    border: 0.01rem solid rgba(226,199,68,1.00);
}



  </style>