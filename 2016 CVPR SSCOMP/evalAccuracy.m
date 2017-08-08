function accuracy = evalAccuracy(label1, label2)
    label2 = bestMap(label1, label2);
    accuracy  = sum(label1(:) == label2(:)) / length(label2);