// src/router/index.js
import { createRouter, createWebHistory } from 'vue-router'
import HomePage from '@/pages/HomePage.vue'
import ContactPage from '@/pages/ContactPage.vue'
import ModuleGeneFinder from '@/pages/ModuleGeneFinder.vue'
import ModuleClinical from '@/pages/ModuleClinical.vue'
import ModuleDepots from '@/pages/ModuleDepots.vue'
import ModuleCellType from '@/pages/ModuleCellType.vue'
import ModuleSingleCell from '@/pages/ModuleSingleCell.vue'
import ModuleSpatial from '@/pages/ModuleSpatial.vue'
import ModulePerturbation from '@/pages/ModulePerturbation.vue'
import ModuleSummary from '@/pages/ModuleSummary.vue'
import VersionTrackerPage from '@/pages/VersionTrackerPage.vue'
import GetStartedPage from '@/pages/GetStartedPage.vue'
import IssuesPage from '@/pages/IssuesPage.vue'

const routes = [
  {
    path: '/',
    name: 'HomePage',
    component: HomePage
  },
  {
    path: '/contact',
    name: 'ContactPage',
    component: ContactPage
  },
  {
    path: '/genefinder',
    name: 'ModuleGeneFinder',
    component: ModuleGeneFinder
  },
  {
    path: '/clinical',
    name: 'ModuleClinical',
    component: ModuleClinical
  } ,
  {
    path: '/depots',
    name: 'ModuleDepots',
    component: ModuleDepots
  },
  {
    path: '/celltype',
    name: 'ModuleCellType',
    component: ModuleCellType
  },
  {
    path: '/singlecell',
    name: 'ModuleSingleCell',
    component: ModuleSingleCell
  },
  {
    path: '/spatial',
    name: 'ModuleSpatial',
    component: ModuleSpatial
  },
  {
    path: '/perturbation',
    name: 'ModulePerturbation',
    component: ModulePerturbation
  },
  {
    path: '/summary',
    name: 'ModuleSummary',
    component: ModuleSummary
  },
  {
    path: '/version',
    name: 'VersionTrackerPage',
    component: VersionTrackerPage
  },
  {
    path: '/started',
    name: 'GetStartedPage',
    component: GetStartedPage
  },
  {
    path: '/issues',
    name: 'IssuesPage',
    component: IssuesPage
  }
]

const router = createRouter({
  history: createWebHistory(),
  routes
})

router.beforeEach((to, from, next) => {
  if (window.gtag) {
    window.gtag('config', 'G-1BCVXDYM6M', {
      page_path: to.path,
    });
  }
  next();
});

export default router
