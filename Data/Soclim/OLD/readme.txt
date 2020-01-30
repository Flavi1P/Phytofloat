Une fois décompressée, le dossier contient des scripts python ainsi que toutes les données nécessaires à la reproduction des figures principales de l'article dans JGR. Dans le dossier '/data' se trouvent les données collectées pendant la campagne (hydro/bio/cyto/comptages avec biovolumes) ainsi que les données des 5 flotteurs (déjà traitées/nettoyées selon la méthode que j'ai préalablement fourni à Stéphane). Deux script sont les plus importants:
- PLS.py  contient la méthode d'apprentissage de la PLS
- prediction.py contient le script d'application de la PLS sur les données flotteurs pour prédire les %Cphyto.

Pour reproduire les figures principales de l'article, il suffit de faire tourner "prediction.py", tous les autres scripts sont appelés automatiquement dedans (Stéphane, tu me confirmera que ça marche ?). Les scripts sont plutôt bien commentés, mais l'ensemble reste complexe. Pour de futures applications, il faudra peut-être repartir de zéro (même si le concept est là).

Les données bateau sont également sous forme de table Excel dans le supplementary data de l'article.

J'espère que ça fonctionne et que ça aidera pour la suite.




import prediction.py
