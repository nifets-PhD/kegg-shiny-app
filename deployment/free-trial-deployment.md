# Free Trial Deployment Options for KEGG Shiny App

## 🆓 Best Free Trial Options (Perfect for 1 Week Testing)

### Option 1: Railway (Recommended for Free Trial) ⭐

**Free Tier**: $5 credit monthly (enough for ~500 hours)
**Trial**: No credit card required for signup
**Perfect for**: 1 week testing

#### Quick Deploy (5 minutes):

1. **Sign up**: https://railway.app (GitHub login, no CC required)
2. **Deploy**: 
   - Click "Deploy from GitHub repo"
   - Connect your `kegg-shiny-app` repository
   - Railway auto-detects Dockerfile
   - Deploy automatically starts

3. **Get URL**: Railway provides a URL like `kegg-shiny-app-production.up.railway.app`

4. **Add to Cloudflare**:
   ```
   Type: CNAME
   Name: kegg
   Target: kegg-shiny-app-production.up.railway.app
   Proxy: Enabled
   ```

**Cost**: FREE for a week! ($5 credit lasts ~3 weeks)

---

### Option 2: Google Cloud Run (Best Free Tier)

**Free Tier**: 2 million requests/month + 400,000 GB-seconds
**Trial**: $300 credit (90 days) - enough for months
**Perfect for**: Long-term free hosting

#### Deploy Steps:

1. **Enable Google Cloud**: https://console.cloud.google.com
2. **Deploy with one command**:
   ```bash
   # Install gcloud CLI first
   gcloud run deploy kegg-shiny-app \
     --source . \
     --platform managed \
     --region us-central1 \
     --allow-unauthenticated \
     --port 3838 \
     --memory 2Gi \
     --cpu 1
   ```

3. **Get URL**: Google provides URL like `kegg-shiny-app-xyz.a.run.app`

**Cost**: FREE indefinitely (within generous limits)

---

### Option 3: Heroku (Classic Choice)

**Free Alternative**: Use Heroku's eco dynos ($5/month, but 1000 free hours for new accounts)
**Trial**: New accounts get credits

#### Deploy Steps:

1. **Install Heroku CLI**: https://devcenter.heroku.com/articles/heroku-cli
2. **Create heroku.yml**:
   ```yaml
   build:
     docker:
       web: Dockerfile
   run:
     web: /usr/bin/shiny-server
   ```

3. **Deploy**:
   ```bash
   heroku create your-kegg-app
   heroku stack:set container
   git push heroku main
   ```

**Cost**: ~$5/month (but new accounts often get free credits)

---

### Option 4: Digital Ocean (3-Month Trial)

**Free Trial**: $200 credit for 60 days (new accounts)
**Perfect for**: Full production testing

#### Steps:
1. Sign up at https://cloud.digitalocean.com
2. Deploy as described in the main guide
3. Use trial credits (lasts months for your app)

**Cost**: FREE with trial credits

---

### Option 5: Render (Great Free Tier)

**Free Tier**: 750 hours/month for web services
**Trial**: No credit card required
**Perfect for**: Development and testing

#### Deploy Steps:

1. **Sign up**: https://render.com (GitHub login)
2. **Create Web Service**:
   - Connect GitHub repo
   - Runtime: Docker
   - Dockerfile path: `./Dockerfile`
   - Port: 3838

3. **Configure**:
   ```yaml
   name: kegg-shiny-app
   env: docker
   plan: free
   dockerfilePath: ./Dockerfile
   envVars:
   - key: PORT
     value: 3838
   ```

**Cost**: FREE (with 750 hours/month limit)

---

## 🚀 Fastest Setup (Under 10 minutes)

### Railway Quick Start:

```bash
# 1. Push your code to GitHub (if not done)
git add .
git commit -m "Ready for Railway deployment"
git push origin main

# 2. Go to https://railway.app
# 3. Click "Deploy from GitHub repo"
# 4. Select kegg-shiny-app repository
# 5. Watch it deploy automatically!
```

### Add Cloudflare (any option):

1. **Add domain to Cloudflare**
2. **Create CNAME record**:
   ```
   Name: kegg
   Target: [your-railway/render/cloudrun-url]
   Proxy: Enabled (orange cloud)
   ```
3. **Enable HTTPS**: Cloudflare handles SSL automatically

---

## 💰 Cost Comparison for 1 Week

| Service | Week 1 Cost | Free Trial | Ease of Setup |
|---------|-------------|------------|---------------|
| Railway | **$0** | $5 credit | ⭐⭐⭐⭐⭐ |
| Google Cloud Run | **$0** | $300 credit | ⭐⭐⭐⭐ |
| Render | **$0** | 750h free | ⭐⭐⭐⭐⭐ |
| Digital Ocean | **$0** | $200 credit | ⭐⭐⭐ |
| Heroku | **$0-5** | Varies | ⭐⭐⭐ |

---

## 🎯 My Recommendation for You

**Use Railway** for your 1-week test:

✅ **Zero setup complexity** - just connect GitHub
✅ **No credit card required** 
✅ **Automatic Docker deployment**
✅ **$5 free credit** (lasts 3+ weeks)
✅ **Easy Cloudflare integration**
✅ **Built-in HTTPS**

### Railway 3-Step Deploy:

1. **Go to**: https://railway.app
2. **Connect**: Your GitHub kegg-shiny-app repo
3. **Done**: Get URL + add to Cloudflare

Takes literally 5 minutes and costs $0 for your week test! 

After your test week, if you want to continue:
- Railway: $5/month
- Google Cloud Run: Likely stays free
- Upgrade to Digital Ocean: $5/month

Want me to walk you through the Railway deployment step-by-step?
