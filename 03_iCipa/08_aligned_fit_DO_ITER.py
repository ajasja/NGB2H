seqs, lm, R2 = rescore_with_new_pos(seqs, lm, fit_fields, weight_string, set_weights, lm_type, extra_cols=extra_cols)

title = str(N_iter)+"_"+fit_type
p = draw_scatter_interactive('Tm', 'score', seqs, y_range=(0, 80),
                                 title=title + '-aligned', save_to_file=False,
                                 tooltips=tooltips + [('pos', '@pos{safe}')])


moved = seqs[seqs.pos != seqs.old_pos]
p.circle('Tm', 'score', fill_alpha=0.8, size=8, source=moved)
bp.show(p)



fit_metric = get_FIT_dataframe(seqs, lm, N_iter, fit_type, N_feat=len(lm.coef_)+1, N_samples=len(seqs))
fit_metrics_all = fit_metrics_all.append(fit_metric)
display(fit_metric)
model_features = get_model_features(lm, fit_fields, N_iter, fit_type)
model_features_all = model_features_all.append(model_features)
display(model_features)


joblib.dump(lm, '{model_out_dir}/{title}.model'.format(**locals()));
df.to_excel('{model_out_dir}/{title}.score.xlsx'.format(**locals()));
bp.save(p, title=title, filename='models\\{title}.plot.html'.format(**locals()), resources=bokeh.resources.INLINE);