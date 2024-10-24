<template>
  <div id="table-container" style="display: none;">
    <table id="geneFinderTable" class="display" style="width:100%; color: black;">
      <thead>
        <tr id="table-headers">
          <!-- Table headers will be dynamically generated based on data -->
        </tr>
      </thead>
      <tbody>
        <!-- Table body will be populated dynamically -->
      </tbody>
    </table>
    <br style="width: 80%;">
    <div style="width: 80%; margin: auto; color: black;">
      <button id="download-data" class="btn btn-primary" style="color: black; max-width: 200px;">Download</button>
    </div>
  </div>
</template>

<script>

</script>

<style scoped>
#geneFinderTable a {
    color: #007bff; /* Bootstrap 默认的超链接颜色 */
    text-decoration: none;
}

/* 只针对 DataTables 表格不换行 */
#geneFinderTable td {
    white-space: nowrap;
}
</style>
