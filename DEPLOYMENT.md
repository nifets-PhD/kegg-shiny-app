# 🚀 Quick Deployment Guide

## Deployment Summary

**Platform**: shinyapps.io (not GitHub Pages - Shiny needs R server)
**Global Cache**: Pre-populated essential pathways 
**Auto-Deploy**: GitHub Actions workflow available

## Step-by-Step Deployment

### 1. Check Requirements
```r
source("deploy.R")
check_deployment_requirements()
```

### 2. Setup shinyapps.io
- Create account: https://www.shinyapps.io/
- Get tokens: https://www.shinyapps.io/admin/#/tokens
- Configure:
```r
rsconnect::setAccountInfo(name='account', token='token', secret='secret')
```

### 3. Deploy with Global Cache
```r
# Deploy with pre-cached essential pathways
deploy_to_shinyapps("kegg-pathway-visualizer", pre_cache = TRUE)
```

Your app will be live at: `https://[account].shinyapps.io/kegg-pathway-visualizer/`

## Why shinyapps.io instead of GitHub Pages?

- **GitHub Pages**: Static files only (HTML/CSS/JS)
- **Shiny Apps**: Need R server execution
- **shinyapps.io**: Official Posit hosting for Shiny apps
- **Free Tier**: 25 active hours/month

## Global Caching Strategy

The app pre-loads these essential pathways for all users:
- Cell cycle, MAPK signaling, PI3K-Akt, AMPK
- Glycolysis, TCA cycle, HIF-1 signaling
- And more common pathways

This ensures fast loading for most users while respecting KEGG servers.
