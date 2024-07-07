<template>
  <div class="app-container">
    <HeaderComponent />
    <main>
      <div class="banner_mod">
        <div class="grid-container_mod">
          <div class="column-1-mod">
            <img src="@/assets/module_images_Gene_finder.png" alt="Gene Finder Module">
          </div>
          <div class="column-2-mod">
            <div class="mod_intro_text">
              <p>This module allows for the identification of genes fulfilling specific criteria. By inputting user-defined thresholds, integrative analyses of the data present in all modules can be investigated.</p>
            </div>
          </div>
        </div>
      </div>
      <hr class="separation_mod">
      
      <div class="button-container">
        <button class="toggle-sidebar-btn" @click="toggleSidebar">
          {{ isSidebarVisible ? 'Hide Sidebar' : 'Show Sidebar' }}
        </button>
      </div>
      <div :class="['container', { 'container-expanded': !isSidebarVisible }]">
        <div v-if="isSidebarVisible" class="sidebar">
          <form @submit.prevent="handleSubmit">
            <div class="form_header">
              Required Input Fields
              <hr class="form_line">
            </div>
            <div class="form_subhead"> > Select P-value cutoff </div>
            <div>
              <input type="radio" id="pValueType_p" name="pValueType" value="p" v-model="vmodel_pValueType" />
              <label for="pValueType_p">unadj. p</label>
              <input type="radio" id="pValueType_FDR" name="pValueType" value="FDR" v-model="vmodel_pValueType" />
              <label for="pValueType_FDR">FDR</label>
            </div>
            <br>
            <label for="pValue">cutoff: </label>
            <input type="text" id="pValue" v-model="vmodel_pValue" class="styled-input">
            
            <div>
              <div class="form_subhead"> > Select cohorts</div>
              <div id="cohort" class="dropdown-check-list" :class="{ visible: isListVisible }" tabindex="100">
                <span class="anchor" @click="toggleListVisibility">Select Cohorts</span>
                <ul class="items">
                  <li v-for="(cohort, index) in cohorts" :key="index">
                    <input type="checkbox" :id="cohort.value" v-model="selectedCohorts" :value="cohort.value" />
                    <label :for="cohort.value">{{ cohort.label }}</label>
                  </li>
                </ul>
              </div>
            </div>
            <div>
              <div v-for="(input, index) in Traits" :key="index">
                <div class="form_subhead"> > Select Trait {{ index + 1 }}</div>
                <select :id="'Trait_' + (index + 1)" v-model="input.selectedOption" class="select-custom">
                  <option value="BMI">BMI</option>
                  <option value="HOMA-IR">HOMA-IR</option>
                  <option value="age">age</option>
                  <option value="WHR">WHR</option>
                  <option value="waist">waist</option>
                  <option value="hip">hip</option>
                  <option value="circ-glucose">circ-glucose</option>
                  <option value="circ-insulin">circ-insulin</option>
                  <option value="circ-TG">circ-TG</option>
                  <option value="circ-cholesterol">circ-cholesterol</option>
                  <option value="circ-HDL">circ-HDL</option>
                  <option value="circ-LDL">circ-LDL</option>
                  <option value="circ-CRP">circ-CRP</option>
                  <option value="HbA1c">HbA1c</option>
                  <option value="WAT LEP secretion">WAT LEP secretion</option>
                  <option value="WAT TNF secretion">WAT TNF secretion</option>
                  <option value="WAT MCP1 secretion">WAT MCP1 secretion</option>
                  <option value="fat cell volume">fat cell volume</option>
                  <option value="basal lipolysis">basal lipolysis</option>
                  <option value="iso lipolysis">iso lipolysis</option>
                  <option value="iso/basal">iso/basal</option>
                </select>
                
                <div class="form_subhead"> > Choose Spearman's Rho Cutoff</div>
                <vue-slider :id="'RangeSlider_' + (index + 1)" 
                             v-model="input.rangeValues" :min="-1" :max="1" :interval="0.01" 
                             style="width: 80%;" :process-style="{ backgroundColor: '#adbec4' }"> 
                </vue-slider>
                
              </div>
              
              <div class="form_subhead"> > Select traits </div>
              <button type="button" @click="addTrait" class="trait-button">+</button> Add Trait
            </div>
            <button type="button" @click="removeLastTrait" class="trait-button">-</button> Remove Last Trait
            <br>
            <br>
            <div class="form_header">
              Additional Input Fields
              <hr class="form_line">
            </div>
            <ToggleButton instanceId="button1" :additionalFieldsConfig="additionalFieldsConfig1" />
            <br>
            <ToggleButton instanceId="button2" :additionalFieldsConfig="additionalFieldsConfig2" />
            <br>
            <div class="form_header">
              Advanced Input Fields
              <hr class="form_line">
            </div>
            <ToggleButton instanceId="button3" :additionalFieldsConfig="additionalFieldsConfig3" />
            <br>
            <div class="form_header">
              Summary of Selected Options
              <hr class="form_line">
            </div>
            <textarea class="text_summary" rows="4" cols="50" v-model="summary"></textarea>
            <br>
            <button type="submit" class="submit_button" id="submit-filter">Submit</button>
            <button id="download-data" class="submit_button" style="display:none;">Download</button>
          </form>
        </div>
        <div class="main-content">
          <div class="tabs">
            <div class="tabs-title">Gene Finder</div>
            <button :class="{ active: currentTab === 'Results' }" @click="setTab('Results')">Results</button>
            <button :class="{ active: currentTab === 'Details' }" @click="setTab('Details')">Details</button>
          </div>
          <div v-if="currentTab === 'Results'">
            <div id="table-container">
              <table id="geneFinderTable" class="display nowrap" style="width:100%;">
                <thead>
                  <tr id="table-headers">
                    <!-- 表头将根据选择的参数动态生成 -->
                  </tr>
                </thead>
                <tbody>
                </tbody>
              </table>
            </div>
          </div>
          <div v-if="currentTab === 'Details'">
            <DetailsComponent :data="detailData" />
          </div>
        </div>
      </div>
    </main>
    <FooterComponent />
  </div>
</template>

<script>
import { ref, computed } from 'vue';
import VueSlider from 'vue-slider-component';
import 'vue-slider-component/theme/material.css';
import ToggleButton from '@/components/ToggleButton.vue';
import DetailsComponent from '@/components/DetailsComponent.vue';
import HeaderComponent from '@/components/HeaderComponent.vue';
import FooterComponent from '@/components/FooterComponent.vue';
import $ from 'jquery';
import 'datatables.net';
import 'datatables.net-buttons';
import 'datatables.net-buttons/js/buttons.html5.js';
import 'jszip';

export default {
  components: {
    HeaderComponent,
    FooterComponent,
    VueSlider,
    ToggleButton,
    DetailsComponent,
  },
  setup() {
    const vmodel_pValue = ref('0.05');
    const vmodel_pValueType = ref('FDR');
    const isListVisible = ref(false);
    const selectedCohorts = ref([]);
    const cohorts = ref([
      { label: "all subcutaneous", value: "all subcutaneous" },
      { label: "all omental", value: "all omental" },
      { label: "Kerr, A. (2020)", value: "Kerr, A. (2020)" },
      { label: "Petrus, P. (2018)", value: "Petrus, P. (2018)" },
      { label: "Arner, P. (2018)", value: "Arner, P. (2018)" },
      { label: "Keller, M. (2017) sc", value: "Keller, M. (2017) sc" },
      { label: "Arner, E. (2012)", value: "Arner, E. (2012)" },
      { label: "Stančáková, A. (2012)", value: "Stančáková, A. (2012)" },
      { label: "Raulerson, C. (2019)", value: "Raulerson, C. (2019)" },
      { label: "Civelek, M. (2017)", value: "Civelek, M. (2017)" },
      { label: "Krieg, L. (2021) sc", value: "Krieg, L. (2021) sc" },
      { label: "Arner, P. (2016) sc", value: "Arner, P. (2016) sc" },
      { label: "Imbert, A. (2022)", value: "Imbert, A. (2022)" },
      { label: "Armenise, C. (2017)", value: "Armenise, C. (2017)" },
      { label: "Winnier, DA. (2015)", value: "Winnier, DA. (2015)" },
      { label: "Nono Nankam, PA. (2020)", value: "Nono Nankam, PA. (2020)" },
      { label: "Vink, RG. (2017)", value: "Vink, RG. (2017)" },
      { label: "Keller, M. (2017) om", value: "Keller, M. (2017) om" },
      { label: "Krieg, L. (2021) om", value: "Krieg, L. (2021) om" },
      { label: "Arner, P. (2016) om", value: "Arner, P. (2016) om" },
      { label: "Barberio, MD. (2019)", value: "Barberio, MD. (2019)" },
      { label: "Proteomics", value: "Proteomics" }
    ]);
    const Traits = ref([]);
    const addTrait = () => Traits.value.push({ selectedOption: '', rangeValues: [-1, 1], value: '' });
    const removeLastTrait = () => {
      if (Traits.value.length > 0) Traits.value.pop();
    };
    const toggleListVisibility = () => isListVisible.value = !isListVisible.value;
    const isSidebarVisible = ref(true);
    const currentTab = ref('Results');
    const resultsData = ref([]);
    const detailData = ref({ description: 'Detailed description here' });
    const additionalFieldsConfig1 = ref([
      {
        id: 'dropdownInputCellType',
        label: 'Cell Type Specificity',
        type: 'dropdown',
        options: [],
        defaultValue: 'no_filter',
      },
    ]);
    const additionalFieldsConfig2 = ref([
      {
        id: 'dropdownInputDepot',
        label: 'Depot Specificity',
        type: 'dropdown',
        options: [
          { label: 'No filtering', value: 'no_filter' },
          { label: 'Subcutaneous', value: 'sc' },
          { label: 'Omental', value: 'om' },
        ],
        defaultValue: 'no_filter',
      },
      {
        id: 'dropdownInputAdipogenesis',
        label: 'Adipogenesis',
        type: 'dropdown',
        options: [
          { label: 'No filtering', value: 'no_filter' },
          { label: 'Not regulated', value: 'Not regulated' },
          { label: 'Down Early', value: 'Down Early' },
          { label: 'Down Intermediate', value: 'Down Intermediate' },
          { label: 'Down Late', value: 'Down Late' },
          { label: 'Transient Up', value: 'Transient Up' },
          { label: 'Transient Down', value: 'Transient Down' },
          { label: 'Transient Mix', value: 'Transient Mix' },
          { label: 'Up Intermediate', value: 'u_inter' },
          { label: 'Up Late', value: 'Up Late' },
        ],
        defaultValue: 'no_filter',
      },
      {
        id: 'dropdownInputSex',
        label: 'Sex Differences',
        type: 'dropdown',
        options: [
          { label: 'No filtering', value: 'no_filter' },
          { label: 'Upregulated in women', value: 'women' },
          { label: 'Upregulated in men', value: 'men' },
        ],
        defaultValue: 'no_filter',
      },
    ]);
    const additionalFieldsConfig3 = ref([
      {
        id: 'dropdownInputGene',
        label: 'Transcript Selection',
        type: 'dropdown',
        options: [
          { label: 'All Genes', value: 'all_genes' },
          { label: 'Only Coding Genes', value: 'protein coding' },
          { label: 'Only Non-Coding Genes', value: 'noncoding' },
        ],
        defaultValue: 'all_genes',
      },
    ]);

    const summary = computed(() => {
      const selectedOptions = [];
      selectedOptions.push(`P-value Cutoff: ${vmodel_pValue.value}, Type: ${vmodel_pValueType.value}`);
      selectedOptions.push(`Selected Cohorts: ${selectedCohorts.value.join(', ')}`);
      Traits.value.forEach((input, index) => {
        selectedOptions.push(`Trait ${index + 1}: ${input.selectedOption}, Rho Cutoff: ${input.rangeValues}`);
      });
      return selectedOptions.join('\n');
    });

    const handleSubmit = () => {
      const requestData = getRequestData();
      initTable(requestData);
      generateCSV(requestData);
    };

    const getRequestData = () => {
      const requestData = {
        Pvaltype: vmodel_pValueType.value,
        PvalCutoff: vmodel_pValue.value,
        cohort: selectedCohorts.value,
      };

      Traits.value.forEach((trait, index) => {
        requestData[`trait${index + 1}`] = trait.selectedOption;
        requestData[`range${index + 1}`] = trait.rangeValues.join(',');
      });

      return requestData;
    };

    const initTable = (requestData) => {
      const ajaxUrl = "https://www.adiposetissue.org/genefinder_DataTables";

      // 检查并销毁现有的DataTable实例
      if ($.fn.DataTable.isDataTable('#geneFinderTable')) {
        $('#geneFinderTable').DataTable().destroy(); // 销毁现有DataTable实例
        $('#geneFinderTable').empty(); // 清空表格内容，包括thead
      }

      $.ajax({
        type: "POST",
        url: ajaxUrl,
        contentType: "application/json",
        data: JSON.stringify(requestData),
        success: (json) => {
          if (json.data.length > 0) {
            const keys = Object.keys(json.data[0]);
            const columns = keys.map(key => ({
              data: key,
              title: key
            }));

            // 重新初始化DataTable
            $('#geneFinderTable').DataTable({
              processing: true,
              serverSide: true,
              ordering: true,
              ajax: {
                type: "POST",
                url: ajaxUrl,
                contentType: "application/json",
                data: (d) => JSON.stringify($.extend({}, d, requestData)),
                dataSrc: "data"
              },
              columns,
              scrollX: true,
              destroy: true,
              initComplete: function () {
                $('#table-container').show();
                this.api().columns.adjust().draw();
              }
            });
          }
        }
      });
    };

    const generateCSV = (requestData) => {
      let csvDownloadUrl;
      let isCsvReady = false;
      $.ajax({
        type: "POST",
        url: "https://www.adiposetissue.org/genefinder_csv",
        contentType: "application/json",
        data: JSON.stringify(requestData),
        success: (response) => {
          csvDownloadUrl = response.download_url;
          isCsvReady = true;
          $('#download-data').show();
          $('#download-data').on('click', () => {
            if (isCsvReady && csvDownloadUrl) {
              window.location.href = csvDownloadUrl;
            } else {
              alert("The CSV file is not ready yet, please try again later.");
            }
          });
        },
        error: (error) => {
          console.log("生成 CSV 文件出错", error);
        }
      });
    };

    const toggleSidebar = () => {
      isSidebarVisible.value = !isSidebarVisible.value;
    };

    const setTab = (tab) => {
      currentTab.value = tab;
    };

    return {
      vmodel_pValue,
      vmodel_pValueType,
      isListVisible,
      selectedCohorts,
      cohorts,
      Traits,
      addTrait,
      removeLastTrait,
      toggleListVisibility,
      isSidebarVisible,
      currentTab,
      resultsData,
      detailData,
      additionalFieldsConfig1,
      additionalFieldsConfig2,
      additionalFieldsConfig3,
      handleSubmit,
      getRequestData,
      initTable,
      generateCSV,
      summary,
      toggleSidebar,
      setTab
    };
  }
};
</script>

<style scoped>
@import url('https://fonts.googleapis.com/css2?family=Red+Hat+Display:wght@400;700&display=swap');

html {
  font-size: 10px; /* 1rem = 10px for easier calculations */
}

body {
  font-family: 'Red Hat Display', sans-serif;
}
h1, h2, h3, h4, h5, h6, button, table {
  font-family: 'Red Hat Display', sans-serif;
}

.container {
  width: 81%;
  height: 100%;
  margin: auto;
  display: flex;
  flex: 1;
  font-family: 'Red Hat Display', sans-serif;
  transition: margin-left 0.3s ease;
}

.container-expanded .main-content {
  margin-left: 0;
}

.sidebar {
  width: 20%;
  padding: 2rem; /* 20px */
  background-color: #EDF2F2;
  border-width:  0.15rem; /* 2px */
  font-family: "Red Hat Display";
  border-color: #1a4659;
  border-radius: 2.5rem; /* 25px */
  border-style: solid;
  overflow-y: auto;
  height: 100%;
  position: sticky;
  top: 0;
  transition: width 0.3s ease, padding 0.3s ease; /* Smooth transition for resizing */
}

.button-container {
  width: 81%;
  margin: 0 auto;
  text-align: left;
}

.toggle-sidebar-btn {
  display: inline-block;
  background-color: #EDF2F2;
  color: #1a4659;
  border-color: #1a4659;
  border-style: solid;
  border-width:  0.15rem; /* 2px */
  padding: 0.6rem 0.7rem;
  cursor: pointer;
  border-radius: 2rem;
  margin-top: 0.5rem;
  margin-bottom: 0.2rem;
  width: 24.5%;
}

.form_header {
  font-weight: bold;
  margin-bottom: 5px;
  margin-top: 10px;
  position: relative;
  display: flex;
  align-items: center;
}
.form_subhead {
  /*text-decoration: underline;*/
  margin-top: 1.2rem; 
  margin-bottom: 0.5rem; /* 5px */
}

/* Hide the default radio button */
input[type="radio"] {
  display: none;
}

/* Style the labels associated with the radio buttons */
input[type="radio"] + label {
  position: relative;
  padding-left: 2rem;
  padding-right: 2rem;
  cursor: pointer;
  line-height: 1rem;
  display: inline-block;
}

/* Create the custom radio button */
input[type="radio"] + label::before {
  content: "";
  position: absolute;
  left: 0;
  top: 0;
  width: 1rem;
  height: 1rem;
  border: 0.08rem solid #b2dfda;
  border-radius: 50%;
  background:  #f7fcfc;
}

input[type="radio"]:checked + label::before {
  background:#1a4659; 
  border-color:#1a4659;
}

/* Create the inner circle for the checked state */
input[type="radio"]:checked + label::after {
  content: "";
  position: absolute;
  left: 5px;
  top: 5px;
  width: 8px;
  height: 8px;
  border-radius: 50%;
  background: #f7fcfc;
}


.styled-input {
  font-family: 'Red Hat Display', sans-serif; 
  border-radius: 1rem; 
  padding: 0.3rem; 
  border: 0.1rem solid #adbec4; 
  background-color:#f7fcfc; 
  font-size: 1rem;
  margin-top: -1rem;
  /*box-shadow: inset 0 1px 3px rgba(0, 0, 0, 0.1);*/
  transition: border-color 0.3s, box-shadow 0.3s;
  width: 58%;
}

.dropdown-check-list {
  display: inline-block;
  font-family: 'Red Hat Display', sans-serif; 
  color:#1a4659;
  font-size: 1rem; 
}

.dropdown-check-list .anchor {
  position: relative;
  cursor: pointer;
  display: inline-block;
  padding: 5px 50px 5px 10px;
  border: 0.1rem solid#adbec4;
  border-radius: 1rem; 
  width: 80%;
}
.dropdown-check-list .anchor:after {
  position: absolute;
  content: "";
  border-left: 0.1rem solid #1a4659;
  border-top: 0.1rem solid #1a4659;
  padding: 0.3rem;
  right: 1rem;
  top: 20%;
  transform: rotate(-135deg);
}
.dropdown-check-list .anchor:active:after {
  right: 1rem;
  top: 20%;
}

.dropdown-check-list ul.items {
  padding: 0.2rem;
  display: none;
  margin: 0;
  border: 0px solid #ccc;
  border-top: none;
}
.dropdown-check-list ul.items li {
  list-style: none;
}
.dropdown-check-list.visible .anchor {
  color:#1a4659;
}
.dropdown-check-list.visible .items {
  display: block;
}

.custom-select option {
  -webkit-appearance: none;
  -moz-appearance: none;
  appearance: none;
}

.select-custom {
  width: 84%; 
  padding: 0.4rem; 
  font-size: 1rem; 
  border: 0.1rem solid #adbec4; 
  border-radius: 2rem; 
  background-color: #f7fcfc; 
  color:#1a4659; 
  appearance: none; 
  cursor: pointer; 
}
/* Style the select dropdown arrow */
.select-custom::after {
  content: '';
  position: absolute;
  right: 10px;
  top: 50%;
  transform: translateY(-50%);
  border-width: 5px;
  border-style: solid;
  border-color: #adbec4 transparent transparent transparent;
  pointer-events: none;
}

/* Customize the dropdown options */
.select-custom option {
  background-color:#f7fcfc;
  color:#1a4659;
  padding: 8px 12px;
}

/* Customize the hover effect for options */
.select-custom option:hover {
  background-color: #1a4659 !important;
  color: #f7fcfc !important;
}

/* Style the scrollbar */
.select-custom::-webkit-scrollbar {
  width: 0.7rem;
  color:  #1a4659;
}

.select-custom::-webkit-scrollbar-track {
  background: #f7fcfc;
}

.select-custom::-webkit-scrollbar-thumb {
  background-color: #f7fcfc;
  border-radius: 6px;
  border: 3px solid #f0f0f0;
}

.select-custom::-webkit-scrollbar-thumb:hover {
  background-color: #1a4659;
}

.text_summary {
  width: 90%;
  font-family: 'Red Hat Display';
  font-size: 1rem; 
  background-color: #f7fcfc;
  border-radius: 1rem;
  border: 0.1rem solid #adbec4;
  padding: 0.5rem;
}

.main-content {
  width: 80%;
  padding-left: 1rem;
  overflow-y: auto;
}
.tabs {
  display: flex;
  align-items: center;
  border-bottom: 0.1rem solid #1a4659;
  margin-bottom: 2rem;
}
.tabs button {
  background: none;
  border: none;
  padding: 1rem 1.5rem;
  margin-right: 1rem;
  cursor: pointer;
  font-size: 1rem;
  color: #1a4659;
}
.tabs button.active {
  background-color: #e0e0e0;
  border-radius: 0rem;
}
.tabs-title {
  font-size: 1rem;
  font-weight: bold;
  color: #1a4659;
  margin-right: 1.6rem;
  /* font-family: 'Red Hat Display', sans-serif; */
}

#submit-filter,
#download-data {
  background-color: #E2C744;
  color: #1a4659;
  border-color: #1a4659;
  padding: 0.7rem 1.4rem;
  cursor: pointer;
  z-index: 1000;
  border-radius: 1rem;
  margin-top: 1rem;
  font-size: 1rem;
}

.trait-button {
  background-color: #1a4659;
  color: #adbec4;
  border: 0.1rem solid #adbec4;
  width: 13.5%;
  height: 10%;
  padding: 0.3rem;
  cursor: pointer;
  border-radius: 50%;
  display: inline-block;
  margin-top: 0.5rem;
  font-weight: 1000;
}

.banner_mod {
  padding-top: 0px;
  vertical-align: middle;
  text-align: center;
  background-color: transparent;
  margin-top: 20px;
}
.banner_mod * {
  align-content: center;
}
.grid-container_mod {
  display: grid;
  vertical-align: middle;
  grid-template-columns: 15% 85%;
  grid-gap: 0rem;
  font-family: "Red Hat Display";
  font-size: 0.875rem;
  color: #1a4659;
  position: relative;
  text-align: center;
}
.column-1-mod, .column-2-mod {
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 1.25rem;
}
.column-1-mod img {
  max-width: 75%;
  max-height: 75%;
  padding-left: 12rem;
}
.mod_intro_text {
  width: 100%;
  padding-left: 2rem;
  padding-right: 8rem;
  padding-top: 0.625rem;
  font-size: 0.875rem;
  color: #1a4659;
  text-align: justify;
}
.separation_mod {
  width: 81%;
  margin-top: rem;
  margin-bottom: 0.1rem;
  border: 0.019rem solid rgba(26, 70, 89, 1.00);
}

</style>
