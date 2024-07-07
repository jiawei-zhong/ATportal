<template>
    <div class="header_col">

<div class="banner">

  <div><router-link to="/"><img src="@/assets/logo_banner_full_2.png" alt="" class="logo_large" /></router-link></div>
    </div>

    <hr class="separation">

    <div class="intro_text"> Open-access to clinical, transcriptional and protein data down to the single-cell level
    </div>
    <form id="myForm" class="search_bar" action="/gene_sum" onsubmit="return modifyAction()" method="get">
    <div class="center-container">
        <div class="search-container">

        <div class="search-bar-container">
            <div class="search-icon">
            <svg xmlns="http://www.w3.org/2000/svg" x="0px" y="0px" width="100" height="100" viewBox="0 0 50 50">
                <path
                d="M 21 3 C 11.621094 3 4 10.621094 4 20 C 4 29.378906 11.621094 37 21 37 C 24.710938 37 28.140625 35.804688 30.9375 33.78125 L 44.09375 46.90625 L 46.90625 44.09375 L 33.90625 31.0625 C 36.460938 28.085938 38 24.222656 38 20 C 38 10.621094 30.378906 3 21 3 Z M 21 5 C 29.296875 5 36 11.703125 36 20 C 36 28.296875 29.296875 35 21 35 C 12.703125 35 6 28.296875 6 20 C 6 11.703125 12.703125 5 21 5 Z">
                </path>
            </svg>
            </div>
            <input type="text" id="gene" name="gene" class="search-bar" placeholder="Search gene, eg. LEP">
        </div>

        <button class="search-button">Search</button>

        </div>
    </div>
    </form>

    </div>
  </template>
  
  <script>
  export default {
    name: 'BannerComponent'
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
    border: px solid #ccc;
    border-radius: 5px;
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