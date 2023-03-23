def error_column(x, l_xerr, r_xerr=None, round=2, sameerr=False):
    col = []
    x = x.astype(float)
    l_xerr = l_xerr.astype(float)
    if(sameerr):
        for _x, _l_xerr in zip(x, l_xerr):
                col.append(f'|{_x:.{round}f}+-<<{_l_xerr:.{round}f}>>|')
    else:
        r_xerr = r_xerr.astype(float)
        for _x, _l_xerr, _r_xerr in zip(x, l_xerr, r_xerr):
            col.append(f'|{_x:.{round}f}++<<+{_r_xerr:.{round}f}>>--<<-{_l_xerr:.{round}f}>>|')
    return col

def latex_formatter(inp):
    return inp.replace('|', '$').replace('+-', '\pm').replace('++', '^').replace('--', '_').replace('<<', '{').replace('>>', '}')

def df_to_latex(df, cols, **kwargs):
    return latex_formatter(df[cols].to_latex(index=False, **kwargs))