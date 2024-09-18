<template>
    <div>
      <h3>Details</h3>
      <p>{{ data.description }}</p>
    </div>
  </template>
  
  <script>
  export default {
    props: {
      data: {
        type: Object,
        required: true,
      },
    },
  };
  </script>
  
  <style scoped>
  h3 {
    margin-bottom: 10px;
  }
  
  p {
    margin-bottom: 10px;
  }
  </style>
  