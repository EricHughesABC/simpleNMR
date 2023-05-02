# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(
    ['simpleNMR.py'],
    pathex=[],
    binaries=[],
    datas=[('cdk-2.7.1.jar', '.'), ('NewTest.class', '.'), ('predictorc.jar', '.'), ('cstables/*', 'csTables'), ('html/results_summary_template.html', 'html'), ('html/w3.css', 'html')],
    hiddenimports=["sklearn.utils._typedefs", "sklearn.neighbors._partition_nodes", "sklearn.utils._weight_vector"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
a.datas += Tree('jre', prefix='jre')
a.datas += Tree('html', prefix='html')
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='simpleNMR',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
