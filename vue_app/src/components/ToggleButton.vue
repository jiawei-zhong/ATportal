<!-- ToggleButton.vue -->
<template>
  <div>
    <button type="button" @click="toggleInputFields" class="button-styles">
      {{ isInputVisible ? '▲ Hide Inputs' : '▼ Show Inputs' }}
    </button>

    <div v-if="isInputVisible">
      <!-- Dynamic radio button for instanceId 'button1' -->
      <div v-if="instanceId === 'button1' && radioValue === 'FACS' || radioValue === 'single cell'">
        <label for="radioInput">Radio Input:</label>
        <input type="radio" id="facs" value="FACS" v-model="radioValue" />
        <label for="facs">FACS</label>
        <input type="radio" id="singleCell" value="single cell" v-model="radioValue" />
        <label for="singleCell">single cell</label>
      </div>

      <!-- Dynamic dropdown menu for instanceId 'button1' -->
      <div v-if="(instanceId === 'button1') && ((radioValue === 'FACS') || (radioValue === 'single cell'))">
        <label for="dropdownInput">Dropdown Input:</label>
        <select :id="additionalFields[0].id" v-model="additionalFields[0].value">
          <option v-for="(option, optionIndex) in getDynamicDropdownOptions()" :key="optionIndex" :value="option.value">
            {{ option.label }}
          </option>
        </select>
      </div>
      <!-- Static dropdowns for instanceId 'button1' and 'button2' -->
      <div v-else>
        <div v-for="(field, index) in additionalFields" :key="index">
          <label :for="field.id">{{ field.label }}:</label>
          <input v-if="field.type !== 'dropdown'" :type="field.type" :id="field.id" v-model="field.value">
          <select v-else-if="field.type === 'dropdown'" :id="field.id" v-model="field.value">
            <option v-for="(option, optionIndex) in field.options" :key="optionIndex" :value="option.value">
              {{ option.label }}
            </option>
          </select>
        </div>
      </div>
    </div>
  </div>
</template>
  
<script>
  import {ref,computed} from 'vue';

  export default {
    props: {
      instanceId: {
        type: String,
        required: true,
      },
      additionalFieldsConfig: {
        type: Array,
        required: true,
      },
    },
    setup(props) {
      const isInputVisible = ref(false);
      const radioValue = ref('FACS');

      const additionalFields = computed(() =>
        props.additionalFieldsConfig.map((additionalField) => ({
          ...additionalField,
          value: additionalField.defaultValue || '',
        }))
      );

      const toggleInputFields = () => {
        isInputVisible.value = !isInputVisible.value;
      };

      // Function to get dynamic dropdown options based on the selected radio button
      const getDynamicDropdownOptions = () => {
        if (radioValue.value === 'FACS') {
          return [
            { label: 'No filtering', value: 'no_filter' },
            { label: 'T, NK & NKT', value: 'FACS_T, NK & NKT' },
            { label: 'B', value: 'FACS_B' },
            { label: 'monocyte & macrophage', value: 'FACS_monocyte & macrophage' },
            { label: 'mast', value: 'FACS_mast' },
            { label: 'FAPs', value: 'FACS_FAPs' },
            { label: 'adipocytes', value: 'FACS_adipocytes' },
            { label: 'vascular', value: 'FACS_vascular' },
          ];
        } else if (radioValue.value === 'single cell') {
          return [
            { label: 'No filtering', value: 'no_filter' },
            { label: 'T, NK & NKT', value: 'Singlecell_T, NK & NKT' },
            { label: 'B', value: 'Singlecell_B' },
            { label: 'monocyte & macrophage', value: 'Singlecell_monocyte & macrophage' },
            { label: 'mast', value: 'Singlecell_mast' },
            { label: 'FAPs', value: 'Singlecell_FAPs' },
            { label: 'adipocytes', value: 'Singlecell_adipocytes' },
            { label: 'vascular', value: 'Singlecell_vascular' },
          ];
        } else {
          // Default options when no radio button is selected or for other radio button values
          return [];
        }
      };

      return {
        isInputVisible,
        radioValue,
        additionalFields,
        toggleInputFields,
        getDynamicDropdownOptions,
      };
    },
  };
</script>

<style scoped>
 .button-styles {
  background-color: #1a4659;
  color: #E2C744;
  border: 1px solid #E2C744;
  padding: 5px 10px;
  cursor: pointer;
  border-radius: 10px;
  display: inline-block; 
  margin-top: 5px; 
}
</style>
  