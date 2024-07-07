<template>
    <header>
      <nav class="navbar_in">
        <router-link to="/">
          <img src="@/assets/logo_navbar_2.png" alt="Adipose Tissue Knowledge Portal Logo" class="logo">
        </router-link>
        <div class="navbar">
            <router-link to="/">Home</router-link>
          <div class="dropdown">
            <button class="dropbtn">Modules <i class="fa fa-caret-down"></i></button>
            <div class="dropdown-content">
                <router-link class="dropdown-item" to="/summary">Summary</router-link>
                <router-link class="dropdown-item" to="/clinical">Clinical Relationships</router-link>
                <router-link class="dropdown-item" to="/depots">Depot Comparisons</router-link>
                <router-link class="dropdown-item" to="/celltype">Cell Type</router-link>
                <router-link class="dropdown-item" to="/singlecell">Single-cell/Single-nucleus RNA Seq.</router-link>
                <router-link class="dropdown-item" to="/spatial">Spatial Transcriptomics</router-link>
                <router-link class="dropdown-item" to="/perturbation">Perturbation</router-link>
                <router-link class="dropdown-item" to="/genefinder">Gene Finder</router-link>
            </div>
          </div>
          <div class="dropdown">
            <button class="dropbtn">About <i class="fa fa-caret-down"></i></button>
            <div class="dropdown-content">
                <router-link class="dropdown-item" to="/started">Get started</router-link>
                <router-link class="dropdown-item" to="/contact">Contribute</router-link>
                <router-link class="dropdown-item" to="/issues">Feedback</router-link>
                <router-link class="dropdown-item" to="/version">Version Updates</router-link>
            </div>
          </div>
        </div>
      </nav>
    </header>
  </template>
  
  <script>
  export default {
    name: 'HeaderComponent',
  };
  </script>
  
  <style scoped>
  header {
    background-color: #1a4659;
  }

  
  .header_col {
        background-color: #1a4659;
        width: auto;
        min-width: 500px;
  }
  
  .logo {
    padding-left: 1.875rem;
    padding-top: 0.625rem;
    padding-bottom: 0.625rem;
    padding-right: 1.18rem;
    width: 35%;
    text-align: left;
  }

  .navbar_in *{
    display: inline-block;
    vertical-align: middle;
    padding-left: auto;
    padding-right: auto;
    margin-left: px;
    margin-right: auto;
}
  
.navbar {
    overflow: hidden;
    background-color: #1a4659;
    font-family: "Red Hat Display";
    float: right;
    padding-top: 0.875rem;
    vertical-align: bottom;
    padding-right: 1.875rem;
}

  /* Links inside the navbar */
.navbar a {
    float: left;
    font-size: 15px;
    color: white;
    text-align: center;
    padding: 0.875rem 1rem;
    text-decoration: none;
}


/* The dropdown container */
.dropdown {
    float: left;
    overflow: hidden;
}

/* Dropdown button */
.dropdown .dropbtn {
    font-size: 0.875rem; 
    border: none;
    outline: none;
    color: white;
    padding: 0.875rem 1rem;
    background-color: inherit;
    font-family: inherit; /* Important for vertical align on mobile phones */
    margin: 0; /* Important for vertical align on mobile phones */
}

/* Add a  background color to navbar links on hover */
.navbar a:hover, .dropdown:hover .dropbtn {
     color: #e2c744;
}

/* Dropdown content (hidden by default) */
.dropdown-content {
    display: none;
    position: absolute;
    background-color: #edf2f2;
    max-width: 10rem;
    box-shadow: 0rem 0.5rem 0.875rem 0rem rgba(0,0,0,0.2);
    z-index: 1;
}

/* Links inside the dropdown */
.dropdown-content a {
    float:left;
    color: #1a4659;
    padding: 0.75rem 1rem;
    display: block;
    text-align: left;
}

/* Add a grey background color to dropdown links on hover */
.dropdown-content a:hover {
     min-width: inherit;
}

/* Show the dropdown menu on hover */
.dropdown:hover .dropdown-content {
     display: block;
} 
  </style>
  