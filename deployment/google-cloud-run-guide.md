# Deploy KEGG Shiny App to Google Cloud Run

## ✅ Why Google Cloud Run is Perfect:
- **60-minute build timeout** (vs Railway's 10 minutes)
- **Free tier**: 2 million requests/month
- **$300 free trial credit** for new accounts
- **Automatic scaling** from 0 to thousands of instances
- **Pay only when used** (likely $0-5/month for your app)

---

## 🚀 Step-by-Step Deployment (15 minutes)

### Step 1: Set Up Google Cloud (5 minutes)

1. **Create Google Cloud Account**:
   - Go to: https://console.cloud.google.com
   - Sign in with Google account
   - Accept $300 free credit (if new account)

2. **Create New Project**:
   - Click "Select Project" → "New Project"
   - Project Name: `kegg-shiny-app`
   - Note the Project ID (e.g., `kegg-shiny-app-123456`)

3. **Enable Required APIs**:
   - Go to "APIs & Services" → "Library"
   - Search and enable:
     - "Cloud Run API"
     - "Cloud Build API" 
     - "Container Registry API"

### Step 2: Install Google Cloud CLI (5 minutes)

Choose your operating system:

**macOS**:
```bash
# Install with Homebrew
brew install --cask google-cloud-sdk

# Or download installer
curl https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-darwin-x86_64.tar.gz | tar -xz
./google-cloud-sdk/install.sh
```

**Windows**:
- Download: https://dl.google.com/dl/cloudsdk/channels/rapid/GoogleCloudSDKInstaller.exe
- Run installer and follow prompts

**Linux**:
```bash
curl https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz | tar -xz
./google-cloud-sdk/install.sh
```

### Step 3: Authenticate and Configure (2 minutes)

```bash
# Authenticate with Google Cloud
gcloud auth login

# Set your project (replace with your actual project ID)
gcloud config set project kegg-shiny-app-123456

# Set default region (London for UK users)
gcloud config set run/region europe-west2

# Verify setup
gcloud config list
```

### Step 4: Deploy Your App (3 minutes)

```bash
# Navigate to your project directory
cd /path/to/kegg-shiny-app

# Deploy directly from source (Cloud Build handles everything)
gcloud run deploy kegg-shiny-app \
  --source . \
  --platform managed \
  --region europe-west2 \
  --allow-unauthenticated \
  --port 3838 \
  --memory 4Gi \
  --cpu 2 \
  --timeout 900 \
  --max-instances 10 \
  --min-instances 0
```

**What this does**:
- Uploads your code to Cloud Build
- Builds Docker image (with 60min timeout)
- Deploys to Cloud Run
- Provides HTTPS URL automatically

---

## 🌐 Step 5: Add to Cloudflare (2 minutes)

After deployment, you'll get a URL like:
`https://kegg-shiny-app-xyz-uc.a.run.app`

1. **Add CNAME Record in Cloudflare**:
   ```
   Type: CNAME
   Name: kegg
   Target: kegg-shiny-app-xyz-uc.a.run.app
   Proxy Status: Enabled (orange cloud)
   ```

2. **SSL/TLS Settings**:
   - Mode: "Full (strict)" 
   - Always Use HTTPS: Enabled

---

## 💰 Cost Breakdown

**Free Tier Includes**:
- 2 million requests/month
- 400,000 GB-seconds/month
- 200,000 CPU-seconds/month

**Your KEGG App Usage** (estimated):
- Light usage: **$0/month** (stays in free tier)
- Medium usage: **$1-5/month**
- Heavy usage: **$5-15/month**

**New Account Bonus**: $300 credit (lasts 6+ months for your app)

---

## 🔧 Monitoring & Management

### View Logs:
```bash
gcloud run services logs tail kegg-shiny-app --region=us-central1
```

### Update App:
```bash
# Just redeploy with same command
gcloud run deploy kegg-shiny-app --source .
```

### Scale Settings:
```bash
# Update instance limits
gcloud run services update kegg-shiny-app \
  --max-instances 20 \
  --min-instances 1 \
  --region us-central1
```

### View Service Info:
```bash
gcloud run services describe kegg-shiny-app --region=us-central1
```

---

## 🎯 Quick Start Commands

Here's the complete sequence once you have gcloud set up:

```bash
# 1. Set project
gcloud config set project YOUR-PROJECT-ID

# 2. Deploy
gcloud run deploy kegg-shiny-app \
  --source . \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 3838 \
  --memory 4Gi \
  --cpu 2

# 3. Get URL
gcloud run services describe kegg-shiny-app \
  --region=us-central1 \
  --format="value(status.url)"
```

---

## 🆘 Troubleshooting

### Build Times Out:
```bash
# Increase Cloud Build timeout
gcloud config set builds/timeout 3600  # 60 minutes
```

### Memory Issues:
```bash
# Increase memory allocation
gcloud run services update kegg-shiny-app \
  --memory 8Gi \
  --region us-central1
```

### Check Build Logs:
```bash
gcloud builds log $(gcloud builds list --limit=1 --format="value(id)")
```

---

## ✅ Advantages Over Railway

- ✅ **60x longer build timeout** (60 min vs 1 min)
- ✅ **True serverless** (scales to zero, pay nothing when idle)
- ✅ **Better performance** (dedicated resources)
- ✅ **Enterprise reliability** (99.95% SLA)
- ✅ **Advanced monitoring** built-in

Ready to deploy? Let me know if you need help with any step!
