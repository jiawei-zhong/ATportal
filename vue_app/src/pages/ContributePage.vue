<template>

  <main>
  <div id="main"> 
      <div id="sidebar_about">
      <div class="sidebar_st">
        <div class="text_sidebar">
              <div><a href="#contribute" class="no-decoration">Contribute to the ATportal</a></div>
              <div><a href="#contact" class="no-decoration">Contact</a></div>
              <div><a href="#cite" class="no-decoration">How to cite the ATportal</a></div>
              <div><a href="#team" class="no-decoration">Team</a></div>
        </div>
      </div>
  </div>
  
    <div id="content" class="vign_articles">
        
            <div class="contact_text" id="contribute">  <h2>Contribute to the ATportal</h2> </div>
        
          <div> <hr class="separation_vign"> </div>
          <div class="contact_text" id="intro_text">
            <div class="intro-container">
              <div class="intro-text">
               <p>
                We invite submissions of datasets, methods, or visualizations relevant to adipose tissue biology
                 to enhance and expand the Adipose Tissue Knowledge Portal. We welcome a broad range of data,
                  including but not limited to clinical datasets with bulk omics (proteomics, transcriptomics),
                   gene expression profiles, and proteomics data characterizing adipose tissue cell types or their
                    perturbations. We also encourage submissions of novel data types not yet included, such as ATAC sequencing or epigenetic studies.
Submissions should include all relevant primary data and statistical analyses, accompanied by a brief description of the dataset
 and its significance to adipose tissue biology, along with references to related publications. Additionally, 
 we welcome the incorporation of novel Shiny modules developed by other labs.
All submissions will be reviewed bi-monthly by the directors and the scientific 
advisory board of the ATportal project. Recommendations for inclusion will be
 followed by a discussion on the specifics of data visualization and integration.
  Please reach out to us via GitHub or the contact form below.
               </p> 
              </div>
              <div class="intro-image">
                <img src="@/assets/images/logo_navbar_2.png" alt="Introduction Image">
              </div>
            </div>
          </div>
  
  <div class="contact_text" id="contact"><h2>Contact</h2>  </div>
    
    <div> <hr class="separation_vign" > </div>
    
   <div class="contact_text">
    <div class="contact-box">
      <form @submit.prevent="submitForm">
      <div class='subheading'>
        Basic Information
        <hr>
      </div>
      
      <div class="form-group">
        <label for="email">Email:</label>
        <input type="email" id="email" v-model="formData.email">
        <span v-if="formErrors.email">{{ formErrors.email }}</span>
      </div>
      
      <div class="form-group">
        <label for="module">Suggested module:</label>
        <input type="text" id="module" v-model="formData.module">
      </div>

      <div class="form-group">
        <label for="clinical">Clinical data available?</label>
        <input type="radio" id="clinical_yes" value="Yes" v-model="formData.clinical">
        <label for="clinical_yes">Yes</label>
        <input type="radio" id="clinical_no" value="No" v-model="formData.clinical">
        <label for="clinical_no">No</label>
      </div>

      <div class="form-group">
        <label for="pmid">PMID:</label>
        <input type="text" id="pmid" v-model="formData.pmid">
      </div>

      <div class='subheading'>
        Data Availability
        <hr>
      </div>

      <div class="form-group">
        <label for="data_type">Data type:</label>
        <input type="radio" id="data_trans" value="Transcriptomics" v-model="formData.dataType">
        <label for="data_trans">Transcriptomics</label>
        <input type="radio" id="data_prot" value="Proteomics" v-model="formData.dataType">
        <label for="data_prot">Proteomics</label>
        <input type="radio" id="data_other" value="Other" v-model="formData.dataType">
        <label for="data_other">Other</label>
      </div>

      <div class="form-group">
        <label for="repository">Repository:</label>
        <input type="text" id="repository" v-model="formData.repository">
      </div>

      <div class="form-group">
        <label for="identifier">Identifier:</label>
        <input type="text" id="identifier" v-model="formData.identifier">
      </div>

      <div class="form-group">
        <label for="message">Message:</label>
        <textarea id="message" v-model="formData.message"></textarea>
      </div>

      <button type="submit">Submit</button>
    </form>
        </div>
   </div>
        
  <div class="contact_text" id="cite"><h2>How to cite the ATportal</h2>  </div>
    
   <div> <hr class="separation_vign" > </div>
    
   <div class="contact_text">
    <p>
      We kindly request that users cite the Adipose Tissue Knowledge Portal,
       as well as the original papers from which the data is derived, when utilizing the content of the portal in their scientific publications.

The ATportal paper is currently available on bioRxiv.
 For individual datasets, please refer to the dataset vignettes under "Get Started." Datasets used in visualizations are also inidcated within the respective modules.
    </p>
   </div>
  
  <div class="contact_text" id="team">  <h2>Team</h2>  </div>
    
    <div> <hr class="separation_vign" > </div>
  
      <div class="contact_text">
        <table class="team-table">
      <tr>
        <td class="team-category">Project directors</td>
        <td>
          <ul class="team-list">
            <li><a href="#">Mikael Rydén</a></li>
            <li><a href="#">Niklas Mejhert</a></li>
            <li><a href="#">Lucas Massier</a></li>
          </ul>
        </td>
      </tr>
      <tr>
        <td class="team-category">Scientific Advisory Board</td>
        <td>
          <ul class="team-list">
            <li>Matthias Blueher, Helmholtz Research Center, Munich, Germany</li>
            <li>Person 2</li>
            <li>Person 3</li>
          </ul>
        </td>
      </tr>
      <tr>
        <td class="team-category">Creative lead</td>
        <td>Danae Zareifi</td>
      </tr>
      <tr>
        <td class="team-category">Software engineers</td>
        <td>
          <ul class="team-list">
            <li>Jiawei Zhong</li>
            <li>Danae Zareifi</li>
            <li>Mattias Hansen</li>
            <li>Lucas Massier</li>
          </ul>
        </td>
      </tr>
      <tr>
        <td class="team-category">Computational biologists</td>
        <td>
          <ul class="team-list">
            <li>Jiawei Zhong</li>
            <li>Danae Zareifi</li>
            <li>Mattias Hansen</li>
            <li>Lucas Massier</li>
            <li>Maheswary Muniandy</li>
          </ul>
        </td>
      </tr>
    </table>
         </div>
    
   <div class="contact_text">
       <p>
      </p>
   </div>
        
    </div>
  </div>
  </main>
  </template>
  <script>
  import emailjs from 'emailjs-com';

  export default {
    name: 'ContributePage',
    data() {
      return {
        formData: {
          email: '',
          module: '',
          clinical: '',
          pmid: '',
          dataType: '',
          repository: '',
          identifier: '',
          message: '',
        },
        formErrors: {}
      };
    },

    methods: {
      submitForm() {
        // Validate email field
        if (!this.formData.email) {
          this.formErrors.email = 'Email is required';
          return;
        }

        const serviceID = 'service_4nr3jof';
        const templateID = 'template_dvhrgpc';
        const userID = 'k8jdo_V1SkI_qpEtp';  

        // Send the form data using emailjs
        emailjs.send(serviceID, templateID, this.formData, userID)
          .then((response) => {
            console.log('Form submitted successfully:', response);
            this.resetForm();  // Reset form on successful submission
          })
          .catch((error) => {
            console.error('Error submitting form:', error);
          });
      },

      resetForm() {
        // Reset form data
        this.formData = {
          email: '',
          module: '',
          clinical: '',
          pmid: '',
          dataType: '',
          repository: '',
          identifier: '',
          message: ''
        };
      }
    }
  }
</script>
 
    <style scoped>
.team-table {
  width: 100%;
  border-collapse: collapse; /* Collapse borders to show lines correctly */
}

.team-category {
  font-weight: bold;
  padding-right: 20px;
  vertical-align: top;
  white-space: nowrap;
}

.team-list {
  list-style: none;
  padding-left: 0;
  margin: 0;
}

.team-list li {
  margin-bottom: 5px;
}

tr {
  border-bottom: 1px solid #5c5858; /* Adds a horizontal line below each row */
}

tr:last-child {
  border-bottom: none; /* Removes the border under the last row */
}

td {
  padding: 15px 0; /* Adds more spacing within each row */
}
    </style>
    