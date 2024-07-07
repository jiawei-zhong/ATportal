<template>
    <footer id="footer">
      <div class="footer_col">
        <div class="grid-container_bot">
          <div class="grid-child">
            <h4>Publication</h4>
            <p>The WATportal manuscript is available via bioRxiv [link to bioRxiv]</p>
          </div>
          <div class="grid-child">
            <h4>Data and Code:</h4>
            <p>The code for this portal is available through <a href="https://github.com/jiawei-zhong/WAT_portal">GitHub</a>.</p>
          </div>
          <div class="grid-child">
            <h4>Get started</h4>
            <p>The <router-link to="/started">Get started</router-link> page provides resources on how to use the portal and contains details on all included data sets</p>
          </div>
          <div class="grid-child">
            <h4>Contribute</h4>
            <p>Check the  <router-link to="/contact">Contribute</router-link> page to learn how to contribute to the WATportal or to contact the team.</p>
          </div>
          <div class="grid-child">
            <h4>Version Update</h4>
            <p><router-link to="/version">{{ latestRelease.name }}</router-link> (updated: {{ formattedDate }})</p>
          </div>
        </div>
      </div>
    </footer>
  </template>
  
  <script>
import axios from 'axios';

export default {
  name: 'FooterComponent',
  data() {
    return {
      latestRelease: {},
    };
  },
  computed: {
    formattedDate() {
      if (this.latestRelease.published_at) {
        const options = { year: 'numeric', month: 'short', day: 'numeric' };
        return new Date(this.latestRelease.published_at).toLocaleDateString('en-US', options);
      }
      return '';
    },
  },
  created() {
    this.fetchLatestRelease();
  },
  methods: {
    async fetchLatestRelease() {
      const owner = 'lmassier';
      const repo = 'hWAT_singlecell';
      try {
        const response = await axios.get(`https://api.github.com/repos/${owner}/${repo}/releases/latest`);
        this.latestRelease = response.data;
      } catch (error) {
        console.error('Error fetching the latest release:', error);
      }
    },
  },
};
  </script>
  
  <style scoped>
  #footer {
        margin-top: auto;
        bottom: 0;
        width: 100%;
        bottom:0;
  }

  .footer_col {
    background-color: #ADBEC4;
    padding-left: 1.875rem;
    padding-top: 3.125rem;
    padding-bottom: 0.938rem;
    margin-top: 0.75rem;
    position: inherit; 
    bottom: 0;
}
  
.grid-container_bot {
    width: 60%;
    display: grid;
    grid-template-columns: 1fr 1fr 1fr 1fr 1fr;
    grid-gap: 0.938rem;
    padding-left: 0.938rem;
    padding-right: 0;
    font-family: "Red Hat Display";
    font-size: 0.75rem;
    color: #1a4659;
    position: relative;
    text-align: left;
}

  .grid-child h4 {
    margin-bottom: 10px;
  }
  
  .grid-child p {
    margin: 0;
  }
  </style>
  