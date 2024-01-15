tld = "./Output/"


def main():
    best_idxs = ["01", "06", "11", "16", "21", "26"]

    folds = os.listdir(raw_data_path)
    for fold in folds:
        with open(tld + fold + "best.dat", "w") as out:
            for idx in best_idxs:
                with open(tld + fold + "best" + idx + ".dat", "r") as inn:
                    lines = inn.readlines()
                    out.writelines(lines)
                    pass
                pass
            pass
        pass

    pass
