import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
font_manager.fontManager.addfont('/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc')
plt.rcParams.update({
    'font.family':'WenQuanYi Zen Hei','font.size':11,'axes.unicode_minus':False,
    'axes.spines.top':False,'axes.spines.right':False,
    'axes.edgecolor':'#4a4a4a','axes.linewidth':0.9,
    'xtick.color':'#4a4a4a','ytick.color':'#4a4a4a','text.color':'#1a1a1a',
    'axes.labelcolor':'#1a1a1a','figure.facecolor':'white','axes.facecolor':'white',
    'savefig.dpi':220,'savefig.bbox':'tight','legend.frameon':False,
})
INK='#1a1a1a'; MUT='#6b6b6b'
BLUE='#2a78d6'; ORANGE='#eb6834'; GREY='#9a9a9a'; TEAL='#1baf7a'; RED='#d03b3b'
# ordinal lesion-stage ramp: control grey, then light->dark blue
STAGE={'control_white_matter':'#9a9a9a','MS_periplaque_white_matter':'#9ec5f4',
       'chronic_inactive_MS_lesion_edge':'#5598e7','chronic_active_MS_lesion_edge':'#2a78d6',
       'MS_lesion_core':'#104281'}
STAGE_CN={'control_white_matter':'对照白质','MS_periplaque_white_matter':'斑块周围白质',
          'chronic_inactive_MS_lesion_edge':'慢性非活动边缘','chronic_active_MS_lesion_edge':'慢性活动边缘',
          'MS_lesion_core':'病灶核心'}
CT_CN={'oligodendrocytes':'少突胶质','astrocytes':'星形胶质','immune':'小胶质/巨噬','neurons':'神经元',
       'opc':'OPC','vascular_cells':'血管/内皮','lymphocytes':'淋巴细胞'}
def violin_strip(ax,groups,colors,labels,ylab,seed=0,maxpts=280,width=0.72):
    rng=np.random.default_rng(seed)
    for i,(v,c) in enumerate(zip(groups,colors)):
        v=np.asarray(v,dtype=float)
        if len(v)>3 and np.ptp(v)>0:
            p=ax.violinplot([v],positions=[i],widths=width,showextrema=False,showmedians=False)
            for b in p['bodies']:
                b.set_facecolor(c); b.set_edgecolor('none'); b.set_alpha(0.32)
        nz=v[v>0]
        if len(nz):
            s=nz if len(nz)<=maxpts else rng.choice(nz,maxpts,replace=False)
            ax.scatter(i+rng.uniform(-0.17,0.17,len(s)),s,s=3.4,color=c,alpha=0.72,linewidths=0)
        if len(v): ax.hlines(np.mean(v),i-0.26,i+0.26,color=c,lw=2.0,zorder=5)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_ylabel(ylab)
