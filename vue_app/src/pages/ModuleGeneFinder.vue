<!-- src/components/HomePage.vue -->
<template>
  <div class="app-container">
    <HeaderComponent />
    <div class="main-container">
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
      <button class="toggle-sidebar-btn" @click="toggleSidebar">
        {{ isSidebarVisible ? 'Hide Sidebar' : 'Show Sidebar' }}
      </button>
      <div :class="['container', { 'container-expanded': !isSidebarVisible }]">
      
      <div v-if="isSidebarVisible" class=sidebar>
        <form>
          <div class="form_header">Required Input Fields <hr class="form_line"> </div>
        
        <!-- p value type -->
          <div class="form_subhead"> P-value cutoffs </div>
          <div>
            <input type="radio" id="pValueType_p" name="pValueType" value="p" v-model="vmodel_pValueType" />
            <label for="pValueType_p">unadj. p</label>
            <input type="radio" id="pValueType_FDR" name="pValueType" value="FDR" v-model="vmodel_pValueType" />
            <label for="pValueType_FDR">FDR</label>
          </div>
        
        <br>

        <!-- p value input -->
        <label for="pValue">cutoff: </label>
        <!-- Use v-model to bind the input value to a data property -->
        <input type="text" id="pValue" v-model="vmodel_pValue">
        <br>

        <!-- Cohorts -->
        <div>
          <div class="form_subhead">Select cohorts</div>
          <div id="cohort" class="dropdown-check-list" :class="{ visible: isListVisible }"  tabindex="100" >
            <span class="anchor" @click="toggleListVisibility">Select Cohorts</span>
            <ul class="items">
              <!-- <li v-for="(cohort, index) in cohorts" :key="index">
                <input type="checkbox" :id="cohort" v-model="selectedCohorts" :value="cohort" />
                <label :for="cohort">{{ cohort }}</label>
              </li> -->
              <li v-for="(cohort, index) in cohorts" :key="index">
                <input type="checkbox" :id="cohort.value" v-model="selectedCohorts" :value="cohort.value" />
                <label label :for="cohort.value">{{ cohort.label }}</label>
              </li>
            </ul>
          </div>
        </div>
        

        <!-- Traits -->
        <div>
          <div v-for="(input, index) in Traits" :key="index">
            <!-- Dropdown menu for the first input -->
          <div class="form_subhead"> Select Trait {{ index + 1 }} </div>
          <select :id="'Trait_' + (index + 1)" v-model="input.selectedOption">
              <option value="BMI">BMI</option> <!-- options should ideally be taken from database-->
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
            <br>
            <!-- Dual range slider -->
            <div class="form_subhead">Choose Spearman's Rho Cutoff</div>
            <vue-slider
              :id="'RangeSlider_' + (index + 1)"
              v-model="input.rangeValues"
              :min="-1"
              :max="1"
              :interval="0.01"
              style="width: 90%;"
              :process-style="{ backgroundColor: '#e8ebe6' }"
              > </vue-slider> 
            <!-- Additional input field -->
            <br>
          </div>
          <br>
          <button type="button" @click="addTrait" class="trait-button">+</button> Add Trait
        </div>
        <button type="button" @click="removeLastTrait" class="trait-button">-</button> Remove Last Trait
        <br>

        
        <div class="form_header">Additional Input Fields
          <hr class="form_line">
        </div>
        <div>
          <!-- Example with different configurations for each instance -->
          <div class="form_subhead">Cell Type Specificity </div>
          <ToggleButton instanceId="button1" :additionalFieldsConfig="[
            // set an empty dropdown, options are rednered in ToggleButton.vue
            { id: 'dropdownInputCellType', 
              label: 'Cell Type Specificity',
              type: 'dropdown',
              options: [],
              defaultValue: 'no_filter',
            },
          ]" />
        </div>
        <br>
        <div>
          <div class="form_subhead">Sex, Depot and Adipogenesis regulation </div>
          <ToggleButton instanceId="button2" :additionalFieldsConfig="[
            { id: 'dropdownInputDepot', 
              label: 'Depot Specificity', 
              type: 'dropdown', 
              options: [
                { label: 'No filtering', value: 'no_filter' },
                { label: 'Subcutaneous', value: 'sc' },
                { label: 'Omental', value: 'om' },
              ],
              defaultValue: 'no_filter',
            },
            { id: 'dropdownInputAdipogenesis', 
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
            { id: 'dropdownInputSex', 
              label: 'Sex Differences', 
              type: 'dropdown', 
              options: [
                { label: 'No filtering', value: 'no_filter' },
                // { label: 'No Sex Specificity', value: 'no_sex' },
                { label: 'Upregulated in women', value: 'women' },
                { label: 'Upregulated in men', value: 'men' },
              ],
              defaultValue: 'no_filter',
            },
          ]" />
          <!-- Add more instances with different configurations as needed -->
        </div>
        <br>
        <div class="form_header">Advanced Input Fields
          <hr class="form_line">
        </div>
        <div>
          <ToggleButton instanceId="button3" :additionalFieldsConfig="[
            { id: 'dropdownInputGene', 
              label: 'Transcript Selection', 
              type: 'dropdown', 
              options: [
                { label: 'All Genes', value: 'all_genes' },
                { label: 'Only Coding Genes', value: 'protein coding' },
                { label: 'Only Non-Coding Genes', value: 'noncoding' },
              ],
              defaultValue: 'all_genes',
            },
            //{ id: 'numericInput', label: 'Numeric Input', type: 'number' },
            //{ id: 'checkboxInput', label: 'Checkbox Input', type: 'checkbox' },
          ]" />
          <!-- Add more instances with different configurations as needed -->
        </div>
        <br>
        <div>
          <!-- Summary of selected options -->
          <div class="form_header">Summary of Selected Options
            <hr class="form_line">
          </div>
          <textarea class="text_summary" rows="4" cols="50" v-model="summary"></textarea>
        </div>

        <br>
        <button class="submit_button" id="submit-filter">Submit</button>
        </form>
      </div>
        <!-- Display the message dynamically 
      <p>Summary of selected options</p>
      <p>Selected Options: {{ selectedOptions }}</p>
      <p>You entered: {{ vmodel_pValue }}</p> -->
      <div class="main-content">
        <!-- Tabs Navigation -->
        <div class="tabs">
          <div class="tabs-title">Gene Finder</div>
          <button :class="{ active: currentTab === 'Results' }" @click="currentTab = 'Results'">Results</button>
          <button :class="{ active: currentTab === 'Details' }" @click="currentTab = 'Details'">Details</button>
        </div>

        <!-- Tabs Content -->
        <div v-if="currentTab === 'Results'">
          <ResultsComponent :data="resultsData" />
        </div>
        <div v-if="currentTab === 'Details'">
          <DetailsComponent :data="detailData" />
        </div>
      </div>
        

      </div>

        <div id="TraitsData" :data-traits="JSON.stringify(Traits)"></div>
      </main>
      </div>
      <FooterComponent />
    </div>
</template>

<script>
//import { createApp, ref } from 'vue';
import { ref } from 'vue';
import VueSlider from 'vue-slider-component';
import 'vue-slider-component/theme/material.css';
import ToggleButton from '@/components/ToggleButton.vue';
import ResultsComponent from '@/components/ResultsComponent.vue';
import DetailsComponent from '@/components/DetailsComponent.vue';
import HeaderComponent from '@/components/HeaderComponent.vue';
import FooterComponent from '@/components/FooterComponent.vue';


export default {
  name: 'ModuleGeneFinder',
  components: {
    HeaderComponent,
    FooterComponent,
    VueSlider,
    ToggleButton,
    ResultsComponent,
    DetailsComponent,
  },
  data() {
    return {
      isSidebarVisible: true,
    };
  },
  computed: {
    containerClass() {
      return {
        'container-expanded': !this.isSidebarVisible,
      };
    },
    summary() {
      const selectedOptions = [];

      // Add selected options from P-value cutoffs
      selectedOptions.push(`P-value Cutoff: ${this.vmodel_pValue}, Type: ${this.vmodel_pValueType}`);

      // Add selected options from Select cohorts
      selectedOptions.push(`Selected Cohorts: ${this.selectedCohorts.join(', ')}`);

      // Add selected options from Traits
      this.Traits.forEach((input, index) => {
        selectedOptions.push(`Trait ${index + 1}: ${input.selectedOption}, Rho Cutoff: ${input.rangeValues}`);
      });

      // Add selected options from additional input fields
      selectedOptions.push(
        ...this.additionalFields
          .map((field) => `${field.label}: ${field.value}`)
          .filter((option) => option !== 'undefined: undefined')
      );
      

      // Combine all selected options into a string
      return selectedOptions.join('\n');
    },
  },
  methods: {
    toggleSidebar() {
      this.isSidebarVisible = !this.isSidebarVisible;
    },
  },
  setup() {
    // Use ref for reactive data
    const vmodel_pValue = ref('0.05');
    const vmodel_pValueType = ref('FDR');
    //const selectedOptions = ref([]);
    const isListVisible = ref(false);
    const selectedCohorts = ref([]);


    const cohorts = [
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
    ];

    const Traits = ref([]);
    const removeLastTrait = () => {
      if (Traits.value.length > 0) {
        Traits.value.pop();
      }
    };
    const toggleListVisibility = () => {
      isListVisible.value = !isListVisible.value;
    };
    const addTrait = () => {
      Traits.value.push({ selectedOption: '', rangeValues: [-1, 1], value: '' });
    }
    const isInputVisible = ref(false);
    const additionalInputValue = ref('');
    const additionalFields = ref([]);
    const toggleTraits = () => {
      isInputVisible.value = !isInputVisible.value;
    };
    const currentTab = ref('Results');
    const resultsData = ref([
      { name: 'Result 1', value: 'Lorem ipsum' },
      { name: 'Result 2', value: 'Dolor sit amet' },
    ]);
    const detailData = ref({
      description: 'Detailed description with lorem ipsum text. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.',
    });
    //const isSidebarVisible = ref(true);
    //const toggleSidebar = () => {
    //  isSidebarVisible.value = !isSidebarVisible.value;
    //};

    // You can directly return the reactive data
    return {
      vmodel_pValue,
      vmodel_pValueType,
      //selectedOptions,
      isListVisible,
      selectedCohorts,
      cohorts,
      additionalFields,
      Traits,
      isInputVisible,
      currentTab,
      resultsData,
      detailData,
      additionalInputValue,
      removeLastTrait,
      toggleTraits,
      addTrait,
      toggleListVisibility,
      //isSidebarVisible,
      //toggleSidebar,
    };
  },
};
</script>

<style scoped>
@import url('https://fonts.googleapis.com/css2?family=Red+Hat+Display:wght@400;700&display=swap');

body {
  font-family: 'Red Hat Display', sans-serif;
}

h1, h2, h3, h4, h5, h6 {
  font-family: 'Red Hat Display', sans-serif;
}

.form_line {
  display: inline-block;
  position: absolute;
  width: 8%;
  margin-left: 10px;
}
.form_header {
  font-weight: bold;
  margin-bottom: 5px;
  margin-top: 10px;
}
.form_subhead {
  text-decoration: underline;
  margin-top: 5px;
  margin-bottom: 5px;
}
.dropdown-check-list {
  display: inline-block;
}

.dropdown-check-list .anchor {
  position: relative;
  cursor: pointer;
  display: inline-block;
  padding: 5px 50px 5px 10px;
  border: 1px solid #ccc;
}

.dropdown-check-list .anchor:after {
  position: absolute;
  content: "";
  border-left: 2px solid black;
  border-top: 2px solid black;
  padding: 5px;
  right: 10px;
  top: 20%;
  -moz-transform: rotate(-135deg);
  -ms-transform: rotate(-135deg);
  -o-transform: rotate(-135deg);
  -webkit-transform: rotate(-135deg);
  transform: rotate(-135deg);
}

.dropdown-check-list .anchor:active:after {
  right: 8px;
  top: 21%;
}

.dropdown-check-list ul.items {
  padding: 2px;
  display: none;
  margin: 0;
  border: 1px solid #ccc;
  border-top: none;
}

.dropdown-check-list ul.items li {
  list-style: none;
}

.dropdown-check-list.visible .anchor {
  color: #0094ff;
}

.dropdown-check-list.visible .items {
  display: block;
}
.text_summary {
  width: 90%;
}


.container {
  display: flex;
  flex: 1;
}
.sidebar {
  width: 20%;
  padding: 20px;
  background-color: #EDF2F2; /* Optional: Background color for the sidebar */
  border-width: 2px;
  font-family: "Red Hat Display";
  border-color: #1a4659;
  border-radius: 25px;
  border-style: solid;
  overflow-y: auto; /* Scroll if the content overflows */
}
.main-content {
  flex: 1; /* Use flex to make it take the remaining space */
  padding: 20px;
  overflow-y: auto; /* Scroll if the content overflows */
}
.tabs {
  display: flex;
  align-items: center;
  border-bottom: 1px solid #1a4659;
  margin-bottom: 20px;
}
.tabs button {
  background: none;
  border: none;
  padding: 10px 20px;
  margin-right: 10px;
  cursor: pointer;
  font-size: 16px;
  color: #1a4659;
}
.tabs button.active {
  font-weight: bold;
  background-color: #e0e0e0;
  border-radius: 5px;
}
.tabs-title {
  font-size: 18px;
  font-weight: bold;
  color: #1a4659;
  margin-right: 20px;
}
.toggle-sidebar-btn {
  display: inline-block; /* Change to inline-block */
  margin-top: 0px; /* Adjust to move the button down */
  margin-bottom: 5px;
  background-color: #E2C744;
  color: #1a4659;
  border-color: #1a4659;
  padding: 10px 20px;
  cursor: pointer;
  border-radius: 5px;
}

#submit-filter {
  background-color: #E2C744;
  color: #1a4659;
  border-color: #1a4659;
  padding: 10px 20px;
  cursor: pointer;
  z-index: 1000;
  border-radius: 5px;
  font-weight: 700;
}
.container-expanded .main-content {
  margin-left: 0;
}

  #button1 {
  background-color: #1a4659;
  color: #E2C744;
  border-color: #E2C744;
  padding: 10px 20px;
  cursor: pointer;
  z-index: 1000;
  border-radius: 5px;
  font-weight: 700;
}

 .trait-button {
  background-color: #1a4659;
  color: #E2C744;
  border: 1px solid #E2C744;
  padding: 5px 5px;
  cursor: pointer;
  border-radius: 10px;
  display: inline-block; 
  margin-top: 5px; 
  font-weight: 900;
}

.dropdown-check-list {
  display: inline-block;
}
.dropdown-check-list .anchor {
  position: relative;
  cursor: pointer;
  display: inline-block;
  padding: 5px 50px 5px 10px;
  border: 1px solid #ccc;
  border-radius: 5px; /* Rounded corners */
  background: white; /* White background */
}
.dropdown-check-list .anchor:after {
  position: absolute;
  content: "";
  border-left: 1px solid black;
  border-top: 1px solid black;
  padding: 5px;
  right: 10px;
  top: 20%;
  -moz-transform: rotate(-135deg);
  -ms-transform: rotate(-135deg);
  -o-transform: rotate(-135deg);
  -webkit-transform: rotate(-135deg);
  transform: rotate(-135deg);
}
.dropdown-check-list .anchor:active:after {
  right: 8px;
  top: 21%;
}
.dropdown-check-list .items {
  z-index: 1;
  position: absolute;
  list-style: none;
  margin: 0;
  padding: 0;
  display: none;
  border: 1px solid #ccc;
  border-top: none;
}
.dropdown-check-list.visible .items {
  display: block;
}
.dropdown-check-list .items li {
  padding: 6px;
  background-color: #fff;
}
.dropdown-check-list .items li:hover {
  background-color: #def;
}

.banner_mod    {
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
    width: 81%;
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
    border: 0.019rem solid rgba(26,70,89,1.00);
}

</style>
