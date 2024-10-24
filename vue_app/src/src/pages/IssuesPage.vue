<template>
    <main>
    <div id="main"> 
    <div id="sidebar_about">
    <div class="sidebar_st">
      <div class="text_sidebar">
        <div><p>Please help us to improve the AT portal by documenting and reporting issues at  <a href="https://github.com/lmassier/hWAT_singlecell/issues">ATportal issues </a> . </p></div>
      </div>
    </div>
</div>

<div id="content" class="vign_articles">

<div class="contact_text"><h2>GitHub Issues</h2>  </div>
  
  <div> <hr class="separation_vign" > </div>
  
 <div class="contact_text">
    <ul v-if="issues.length">
        <li v-for="issue in issues" :key="issue.id">
          <h2>{{ issue.title }}</h2>
          <p>Created by: {{ issue.user.login }} on {{ new Date(issue.created_at).toLocaleDateString() }}</p>
          <p>{{ issue.body }}</p>
          <a :href="issue.html_url" target="_blank">View Issue</a>
        </li>
      </ul>
      <p v-else>No issues found.</p>
 </div>
</div>
</div>
</main>
  </template>
  
  <script>
  import axios from 'axios';
  
  export default {
    name: 'IssuesPage',
    data() {
      return {
        issues: [],
      };
    },
    created() {
      this.fetchIssues();
    },
    methods: {
      async fetchIssues() {
        const owner = 'jiawei-zhong';
        const repo = 'ATportal';
        try {
          const response = await axios.get(`https://api.github.com/repos/${owner}/${repo}/issues`);
          this.issues = response.data;
        } catch (error) {
          console.error('Error fetching issues:', error);
        }
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
    /*display: block;*/
    padding: 0.5rem 0;
    white-space: nowrap;
    text-decoration: underline;
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
  