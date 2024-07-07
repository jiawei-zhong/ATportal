// src/main.js
import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import HeaderComponent from '@/components/HeaderComponent.vue'
import FooterComponent from '@/components/FooterComponent.vue'
import './assets/styles/global.css'  // Import the global stylesheet

const app = createApp(App)

app.component('HeaderComponent', HeaderComponent)
app.component('FooterComponent', FooterComponent)

app.use(router)
app.mount('#app')