
        .compile -v idl//dotf.pro
        .run idl//dotf.pro
        dotf, 'idl/time.txt', 'idl/flux.txt','1','idl/tf_lnp.sav'
        exit
        