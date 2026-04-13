#!/usr/bin/env python3
"""
ppigFinder — Protein-Protein Interaction Genomic Finder
Entry point for the application.
"""
import sys
import os
import logging
import argparse

# --- O CONSERTO DO CAMINHO ---
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

def setup_logging(debug_mode):
    """Configura o nível de fofoca do terminal."""
    level = logging.DEBUG if debug_mode else logging.INFO
    # Agora o log vai mostrar: [Hora] - [Nome do Arquivo] - [Nível]: Mensagem
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    )

def main():
    # Adiciona o suporte a argumentos de terminal
    parser = argparse.ArgumentParser(description="ppigFinder GUI Application")
    parser.add_argument('--debug', action='store_true', help='Ativa logs detalhados no terminal')
    # Permite passar argumentos extras que o PyQt possa precisar
    args, unknown = parser.parse_known_args()

    setup_logging(args.debug)
    logger = logging.getLogger(__name__)
    logger.info("Starting ppigFinder...")
    
    try:
        from ppigfinder.gui.app import QApplication, ppigFinderApp, _setup_emoji_font, QT_VERSION
        
        # Passa o sys.argv (incluindo os argumentos desconhecidos) pro Qt
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
