function v_chopped = threshold(v, cutoff)

v_chopped = v;
v_chopped(v > mean(v) + cutoff*std(v)) = mean(v) + cutoff*std(v);
v_chopped(v < mean(v) - cutoff*std(v)) = mean(v) - cutoff*std(v);