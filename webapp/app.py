"""
OpenFront Map Generator - Flask Backend
========================================
A web API for generating styled terrain maps from DEM data.
"""

import os
import json
import tempfile
import shutil
from datetime import datetime
from flask import Flask, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
from dotenv import load_dotenv
from functools import wraps

# Load environment variables
load_dotenv()

# Import our map generator
from map_processor import MapProcessor

app = Flask(__name__, static_folder='static')
CORS(app)  # Allow cross-origin requests from GitHub Pages

# Configuration
OPENTOPO_API_KEY = os.environ.get('OPENTOPO_API_KEY', '')
SUPABASE_URL = os.environ.get('SUPABASE_URL', '')
SUPABASE_ANON_KEY = os.environ.get('SUPABASE_ANON_KEY', '')
SUPABASE_SERVICE_KEY = os.environ.get('SUPABASE_SERVICE_KEY', '')

# Output directory for generated maps
OUTPUT_DIR = os.path.join(tempfile.gettempdir(), 'openfront_maps')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def require_auth(f):
    """Decorator to require authentication."""
    @wraps(f)
    def decorated(*args, **kwargs):
        auth_header = request.headers.get('Authorization', '')
        if not auth_header.startswith('Bearer '):
            return jsonify({'error': 'Missing or invalid authorization header'}), 401
        
        token = auth_header.split(' ')[1]
        
        # For this app, we trust the frontend auth (Supabase handles it client-side)
        # The real protection is that users must provide their own API key
        # If you need server-side verification, use Supabase JWT verification
        if token and len(token) > 10:
            request.user = {'token': token}
        elif not SUPABASE_URL:
            # Dev mode - no auth required
            request.user = {'id': 'dev-user'}
        else:
            return jsonify({'error': 'Invalid token'}), 401
        
        return f(*args, **kwargs)
    return decorated


@app.route('/')
def index():
    """Serve the main application page."""
    return send_from_directory('static', 'index.html')


@app.route('/api/health')
def health():
    """Health check endpoint for Render."""
    return jsonify({
        'status': 'healthy',
        'api_key_configured': bool(OPENTOPO_API_KEY),
        'supabase_configured': bool(SUPABASE_URL)
    })


@app.route('/api/generate', methods=['POST'])
@require_auth
def generate_map():
    """
    Generate a styled terrain map.
    
    Request body:
    {
        "name": "My Map",
        "bounds": {
            "south": 35.0,
            "west": 32.0,
            "north": 36.0,
            "east": 34.0
        },
        "dem_source": "COP90"
    }
    """
    try:
        data = request.get_json()
        
        # Validate input
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        name = data.get('name', f'map_{datetime.now().strftime("%Y%m%d_%H%M%S")}')
        bounds = data.get('bounds')
        dem_source = data.get('dem_source', 'COP90')
        width_px = data.get('width_px')
        height_px = data.get('height_px')
        
        if not bounds:
            return jsonify({'error': 'Bounds are required'}), 400
        
        required_bounds = ['south', 'west', 'north', 'east']
        for key in required_bounds:
            if key not in bounds:
                return jsonify({'error': f'Missing bound: {key}'}), 400
        
        # Validate bounds
        south, west, north, east = bounds['south'], bounds['west'], bounds['north'], bounds['east']
        if south >= north or west >= east:
            return jsonify({'error': 'Invalid bounds'}), 400
        
        # Get API key - prefer user-provided, fallback to server config
        user_api_key = data.get('api_key', '')
        api_key = user_api_key if user_api_key else OPENTOPO_API_KEY
        
        if not api_key:
            return jsonify({'error': 'No API key provided. Please set your OpenTopography API key in Settings.'}), 400
        
        # Create output directory for this map
        map_id = f"{name.lower().replace(' ', '_')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        map_dir = os.path.join(OUTPUT_DIR, map_id)
        os.makedirs(map_dir, exist_ok=True)
        
        # Process the map
        processor = MapProcessor(
            api_key=api_key,
            output_dir=map_dir
        )
        
        result = processor.generate(
            name=name,
            south=south,
            west=west,
            north=north,
            east=east,
            width_px=width_px,
            height_px=height_px,
            dem_source=dem_source
        )
        
        return jsonify({
            'success': True,
            'map_id': map_id,
            'files': result['files'],
            'download_url': f'/api/download/{map_id}'
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/download/<map_id>')
def download_map(map_id):
    """Download generated map files as a ZIP."""
    # Note: Auth removed for direct browser downloads
    # In production, you may want to use signed URLs or session-based auth
    map_dir = os.path.join(OUTPUT_DIR, map_id)
    
    if not os.path.exists(map_dir):
        return jsonify({'error': 'Map not found'}), 404
    
    # Create a ZIP file
    zip_path = os.path.join(OUTPUT_DIR, f'{map_id}.zip')
    shutil.make_archive(zip_path.replace('.zip', ''), 'zip', map_dir)
    
    return send_file(zip_path, as_attachment=True, download_name=f'{map_id}.zip')


@app.route('/api/download/<map_id>/<filename>')
def download_file(map_id, filename):
    """Download a specific file from a generated map."""
    map_dir = os.path.join(OUTPUT_DIR, map_id)
    
    if not os.path.exists(map_dir):
        return jsonify({'error': 'Map not found'}), 404
    
    file_path = os.path.join(map_dir, filename)
    if not os.path.exists(file_path):
        return jsonify({'error': 'File not found'}), 404
    
    return send_file(file_path, as_attachment=True)


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)
