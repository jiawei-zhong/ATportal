<template>
  <div>
    <div class="app-container">
      <main>
        <div class="iframe-container">
          <!-- Add a 'ref' attribute to the iframe for easy reference in Vue -->
          <iframe ref="myIframe" height="800" frameborder="0" src="https://atportal.medh.ki.se/shiny/Module_Singlecell/"></iframe>
        </div>
      </main>
    </div>
  </div>
</template>

<script>
import { iframeResizer } from 'iframe-resizer'; // Ensure correct import of iframeResizer

export default {
  name: 'ModuleSingleCell',
  mounted() {
    // Use this.$refs to access the iframe and initialize iframe-resizer
    this.$nextTick(() => {
      iframeResizer({
        log: true, // Enable logging for debugging
        checkOrigin: false // Set to false if iframe content is from the same origin
      }, this.$refs.myIframe);
    });
  }
}
</script>

<style scoped>

</style>
  
  