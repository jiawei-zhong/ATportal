<template>
    <main>
<div id="main"> 
    <div id="sidebar_about">
    <div class="sidebar_st">
      <div class="text_sidebar">
        <div v-for="release in releases" :key="release.id">
              <a :href="`#release-${release.id}`">{{ release.name }}</a>
            </div>
      </div>
    </div>
</div>

<div id="content" class="vign_articles">

<div class="contact_text"><h2>GitHub Releases</h2>  </div>
  
  <div> <hr class="separation_vign" > </div>
  
 <div class="contact_text">
      <ul v-if="releases.length">
        <li v-for="release in releases" :key="release.id">
          <h2>{{ release.name }}</h2>
          <p>{{ release.body }}</p>
          <p>Published at: {{ release.published_at }}</p>
          <a :href="release.html_url" target="_blank">View Release</a>
        </li>
      </ul>
      <p v-else>No releases found.</p>
 </div>
</div>
</div>
</main>

    
    
    
    <div>
    

    </div>
  </template>
  
  <script>
  import axios from 'axios';
  
  export default {
  name: 'VersionTrackerPage',
  data() {
    return {
      releases: [],
      selectedRelease: null,
    };
  },
  created() {
    this.fetchReleases();
  },
  methods: {
    async fetchReleases() {
      const owner = 'jiawei-zhong';
      const repo = 'ATportal';
      try {
        const response = await axios.get(`https://api.github.com/repos/${owner}/${repo}/releases`);
        this.releases = response.data;
        if (this.releases.length > 0) {
          this.selectedRelease = this.releases[0];
        }
      } catch (error) {
        console.error('Error fetching releases:', error);
      }
    },
    showRelease(release) {
      this.selectedRelease = release;
    },
  },
};
  </script>
  
  <style scoped>
  #sidebar {
    display: table-cell;
    width: 25%;
    vertical-align: top;
    color: #adbec4;
    padding-top: 1.25rem;
  }
  
  #content {
    display: table-cell;
    width: 75%;
    vertical-align: top;
  }
  
  .sidebar_st {
    border-radius: 1.563rem;
    background: #edf2f2;
    padding: 0;
    width: 100%;
    height: auto;
  }
  
  .text_sidebar {
    text-align: justify;
    padding-left: 1.25rem;
    padding-right: 1.25rem;
    padding-top: 0.938rem;
    padding-bottom: 1.25rem;
    color: #1a4659;
    font-family: "Red Hat Display";
    font-size: 1rem;
  }
  
  .text_sidebar a {
    color: #1a4659;
    text-decoration: none;
    display: block;
    padding: 0.5rem 0;
  }
  
  .contact_text {
    width: 100%;
    padding-left: 40px;
    padding-top: 0px;
    font-size: 14px;
    color: #1a4659;
    text-align: justify;
    vertical-align: top;
  }
  
  .separation_cont {
    max-width: 200px;
    float: left;
    border: 0.5px solid rgba(26,70,89,1.00);
  }
  
  #main {
    padding-top: 30px;
    padding-left: 20px;
    padding-right: 20px;
    text-align: center;
  }
  
  .sidebar {
    position: static;
    height: auto;
    width: 150px;
    top: 0px;
    float: right;
    margin-top: 100px;
    padding-top: 40px;
    padding-left: 20px;
    background-color: lightblue;
  }
  
  .sidebar div {
    padding: 8px;
    font-size: 24px;
    display: block;
  }
  
  .vign_articles {
    margin-right: 150px;
    max-width: 70%;
    padding-left: 20px;
    padding-top: 20px;
    padding-bottom: 20px;
    padding-right: 20px;
    margin-top: 0px;
    color: #1a4659;
  }
  
  .separation_vign {
    max-width: 1210px;
    float: inherit;
    border: 1px solid rgba(226,198,68,1.00);
    margin-left: 40px;
    margin-top: -15px;
    border-radius: 2px;
  }
  </style>
  