# OpenFront Map Generator - Web App

A polished web application for generating styled terrain maps from DEM data.

## Features

- 🗺️ Interactive map selection with Leaflet
- 🎨 OpenFront-styled terrain color palette
- 🏔️ Automatic DEM downloading from OpenTopography
- 🌊 Rivers and lakes overlay from Natural Earth
- 📍 Automatic province/country detection
- 👤 User authentication with Supabase
- ☁️ Deployable to Render (free tier)

## Quick Start (Local Development)

### 1. Install Dependencies

```bash
cd webapp
pip install -r requirements.txt
```

### 2. Set Up Environment Variables

Copy the example env file and add your API key:

```bash
cp .env.example .env
```

Edit `.env` and add your OpenTopography API key:
```
OPENTOPO_API_KEY=your_key_here
```

Get a free API key at: https://portal.opentopography.org/myopentopo

### 3. Run Locally

```bash
python app.py
```

Open http://localhost:5000 in your browser.

**Note:** In local development mode, authentication is skipped.

---

## Deployment to Render (Free Tier)

### Step 1: Create a Supabase Project (Free)

1. Go to https://supabase.com and sign up
2. Create a new project
3. Go to **Settings → API**
4. Copy these values:
   - `Project URL` → SUPABASE_URL
   - `anon public` key → SUPABASE_ANON_KEY
   - `service_role` key → SUPABASE_SERVICE_KEY

### Step 2: Update Frontend Config

Edit `webapp/static/index.html` and update the CONFIG section:

```javascript
const CONFIG = {
    API_URL: 'https://your-app-name.onrender.com',
    SUPABASE_URL: 'https://your-project.supabase.co',
    SUPABASE_ANON_KEY: 'your-anon-key'
};
```

### Step 3: Push to GitHub

```bash
git add webapp/
git commit -m "Add web application"
git push
```

### Step 4: Deploy to Render

1. Go to https://render.com and sign up
2. Click **New → Web Service**
3. Connect your GitHub repository
4. Configure:
   - **Name:** openfront-map-generator
   - **Root Directory:** webapp
   - **Runtime:** Python 3
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `gunicorn app:app`

5. Add Environment Variables:
   - `OPENTOPO_API_KEY` - Your OpenTopography API key
   - `SUPABASE_URL` - Your Supabase project URL
   - `SUPABASE_ANON_KEY` - Your Supabase anon key
   - `SUPABASE_SERVICE_KEY` - Your Supabase service key
   - `FLASK_SECRET_KEY` - A random string for session security

6. Click **Create Web Service**

### Step 5: Host Frontend on GitHub Pages (Optional)

For faster loading, you can host the frontend separately:

1. Copy `webapp/static/index.html` to a new repo
2. Enable GitHub Pages in repo settings
3. Update `API_URL` in the HTML to point to your Render backend

---

## Project Structure

```
webapp/
├── app.py              # Flask API server
├── map_processor.py    # DEM processing (replaces QGIS)
├── requirements.txt    # Python dependencies
├── render.yaml         # Render deployment config
├── .env.example        # Environment variables template
└── static/
    └── index.html      # Frontend (single-page app)
```

---

## API Endpoints

### `GET /api/health`
Health check endpoint

### `POST /api/generate`
Generate a new map

**Request:**
```json
{
    "name": "Cyprus",
    "bounds": {
        "south": 34.5,
        "west": 32.0,
        "north": 35.7,
        "east": 34.6
    },
    "width": 2048,
    "dem_source": "COP90"
}
```

**Response:**
```json
{
    "success": true,
    "map_id": "cyprus_20260109_143022",
    "files": ["cyprus.png", "cyprus.json"],
    "download_url": "/api/download/cyprus_20260109_143022"
}
```

### `GET /api/download/<map_id>`
Download generated map as ZIP

---

## Limitations

- Free Render tier spins down after 15 minutes of inactivity
- First request after spin-down takes ~30 seconds
- Large maps (4096px+) may timeout on free tier
- OpenTopography API has rate limits

---

## Troubleshooting

### "Host requires authentication" error
→ Make sure your `OPENTOPO_API_KEY` is set correctly

### Map generation times out
→ Try a smaller area or lower resolution (1024px)

### Authentication issues
→ Check your Supabase keys are correctly configured

---

## License

MIT License - Use freely for your projects!
