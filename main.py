#!/usr/bin/env python3
"""
ppigFinder — Protein-Protein Interaction Genomic Finder
Entry point for the application.
"""
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def main():
    logger = logging.getLogger(__name__)
    logger.info("Starting ppigFinder...")
    
    try:
        # Importamos as ferramentas gráficas da nossa nova gaveta
        from ppigfinder.gui.app import QApplication, ppigFinderApp, _setup_emoji_font, QT_VERSION
        
        app = QApplication(sys.argv)
        app.setApplicationName('ppigFinder')
        app.setApplicationDisplayName('ppigFinder — Protein-Protein Interaction Genomic Finder')
        app.setApplicationVersion('1.01')
        app.setStyle('Fusion')
        
        _setup_emoji_font(app)
        
        window = ppigFinderApp()
        window.show()
        
        logger.info("GUI loaded successfully. Running event loop...")
        sys.exit(app.exec() if QT_VERSION == 6 else app.exec_())
        
    except ImportError as e:
        logger.error(f"Failed to load the GUI: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
