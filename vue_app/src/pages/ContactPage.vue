<template>
  <div>
    <div class="app-container">
      <!--
      <div class="contact-page">
        <div class="contact-box">
          <h1>Contact/ Feedback</h1>
          
          <form @submit.prevent="submitForm">
            <div class="form-group">
              <label for="name1">Name:</label>
              <input type="text" id="name1" v-model="form.name" required>
            </div>
            <div class="form-group">
              <label for="email1">Email:</label>
              <input type="email" id="email1" v-model="form.email" required>
            </div>
            <div class="form-group">
              <label for="message1">Message:</label>
              <textarea id="message1" v-model="form.message" required></textarea>
            </div>
            <button type="submit">Submit</button>
          </form>
          <p>{{ formStatus }}</p>
        </div>

        <div class="contact-box">
          <h1>Contribute New Data </h1>
          
          <form>
            <div class='subheading'>
                Basic Information
                <hr>
              </div>
            <div class="form-group">
              <label for="module">Suggested module:</label>
              <input type="text" id="module" name="module">
            </div>
            <div class="form-group">
              <label for="clinical">Clinical data available?</label>
              <input type="radio" id="clinical_yes" name="clinical" value="clinical_yes">
                <label for="clincal_yes">Yes</label>
              <input type="radio" id="clinical_no" name="clinical" value="clinical_no">
                <label for="clincal_no">No</label>
            </div>
            <div class="form-group">
              <label for="pmid">PMID:</label>
              <input type="text" id="pmid" name="pmid">
            </div>
            <div class='subheading'>
                Data availability
                <hr>
              </div>
              <div class="form-group">
              <label for="data_type">Data type</label>
              <input type="radio" id="data_trans" name="data_type" value="data_trans">
                <label for="data_trans">Transcriptomics</label>
              <input type="radio" id="data_prot" name="data_type" value="data_prot">
                <label for="data_prot">Proteomics</label>
              <input type="radio" id="data_other" name="data_type" value="data_other">
                <label for="data_other">Other</label>
            </div>
            <div class="form-group">
              <label for="repository">Repository</label>
              <input type="text" id="repository" name="repository">
            </div>
            <div class="form-group">
              <label for="identifier">Identifier</label>
              <input type="text" id="identifier" name="identifier">
            </div>
            <div class="form-group">
              <label for="message2">Message:</label>
              <textarea id="message2" name="message2"></textarea>
            </div>
            <button type="submit">Submit</button>
          </form>
        </div>
      </div>
    -->
    </div>
  </div>
</template>

<script>
export default {
  name: 'ContactPage',
  data() {
    return {
      form: {
        name: '',
        email: '',
        message: ''
      },
      formStatus: ''
    };
  },
  methods: {
    async submitForm() {
      try {
        const response = await fetch('/api/contact', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify(this.form)
        });

        if (response.ok) {
          this.formStatus = 'Message sent successfully!';
          this.form.name = '';
          this.form.email = '';
          this.form.message = '';
        } else {
          this.formStatus = 'Failed to send message.';
        }
      } catch (error) {
        this.formStatus = 'Error sending message.';
      }
    }
  }
}
</script>

<style scoped>
.app-container {
  display: flex;
  flex-direction: column;
  align-items: center;
  max-width: 1600px; /* Adjust max width as needed */
  margin: 0 auto; /* Center the content horizontally */
}

.contact-page {
  display: flex;
  gap: 30px; /* Adjust the space between the boxes */
  justify-content: space-between; /* Distribute space between boxes */
  margin-top: 20px; /* Add margin to separate from header */
}

.contact-box {
  flex: 1; /* Each box takes equal space */
  padding: 20px;
  background-color: #EDF2F2; /* Background color for the boxes */
  border-width: 2px;
  font-family: "Red Hat Display";
  border-color: #1a4659;
  border-radius: 25px;
  border-style: solid;
  overflow-y: auto; /* Scroll if the content overflows */
  box-sizing: border-box; /* Include padding and border in the element's total width and height */
}

.contact-box h1 {
  text-align: center; /* Center the title */
  margin-bottom: 20px; /* Add space below the title */
}

.contact-box form {
  margin-bottom: 20px; /* Add space below the form */
}

.contact-box p {
  margin-top: 20px; /* Add space above the content paragraphs */
}

/* Form Group Styling */
.form-group {
  display: flex;
  align-items: center;
  margin-bottom: 15px; /* Add space between form groups */
}

.form-group label {
  width: 100px; /* Fixed width for labels */
  margin-right: 10px; /* Space between label and input */
  text-align: left; /* Align label text to the left */
}

.form-group input,
.form-group textarea {
  flex: 1; /* Input fields take the remaining space */
  padding: 10px;
  border: 1px solid #ccc;
  border-radius: 5px;
  font-family: "Red Hat Display";
  text-align: right; /* Align input text to the right */
}

/* Input and Textarea Styling */
input[type="text"],
input[type="email"],
textarea {
  width: 100%; /* Full width within the flex container */
}

textarea {
  height: 100px; /* Fixed height for textarea */
  resize: vertical; /* Allow vertical resizing */
}

/* Button Styling */
button {
  padding: 10px 20px;
  border: none;
  border-radius: 5px;
  background-color: #1a4659;
  color: white;
  cursor: pointer;
  font-family: "Red Hat Display";
}

button:hover {
  background-color: #123244; /* Darker shade on hover */
}
</style>
