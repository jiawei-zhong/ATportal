<template>
    <div>
      <div class="app-container">
        <BannerComponent />
        <CircleComponent />
      </div>
    </div>
  </template>

<script>
import BannerComponent from '@/components/BannerComponent.vue';
import CircleComponent from '@/components/CircleComponent.vue';

export default {
  name: 'HomePage',
  components: {
    BannerComponent,
    CircleComponent
  }
}
</script>