import re


def parser_for_output_file_name(complete_file_name, typ):
    # extract relevant fields from complete_file_name
    map__field__value = {}
    map__field__value["complete_file_name"] = complete_file_name

    first_index = 1 + complete_file_name.rfind("/")
    first_index = 0 if first_index == 0 else first_index
    new_complete_file_name = complete_file_name[first_index:]

    if typ.lower() == 'nudhy':
        new_complete_file_name = new_complete_file_name.replace("__", "/")
        tokens = new_complete_file_name.split("/")
        dataset = tokens[0]
        map__field__value["dataset"] = dataset
        for index in range(1, len(tokens) - 1):
            c_token = tokens[index]
            separator_index = 1 + c_token.find("_")
            c_field_name = c_token[:separator_index - 1]
            c_field_value = c_token[separator_index:]
            map__field__value[c_field_name] = c_field_value
            if c_field_name == "bias":
                c_field_value = float(c_field_value)
                map__field__value[c_field_name] = c_field_value
                continue
            if c_field_name != "algorithm":
                c_field_value = int(c_field_value)
                map__field__value[c_field_name] = c_field_value
                continue
            map__field__value[c_field_name] = c_field_value
    else:
        # CASE: Base, BaseD, ReDi, ReDiD, UnpretentiousNullModel
        tokens = new_complete_file_name.split("_" + typ)
        map__field__value["dataset"] = tokens[0]
        tokens2 = tokens[1].split("_")
        map__field__value["randomSeed"] = extract_seed(tokens2[0])
        map__field__value["algorithm"] = typ
    return map__field__value


def extract_seed(input_string):
    # Extract the seed from input_string
    pattern = r'\d+'  # one or more digits

    matches = re.findall(pattern, input_string)
    if matches:
        return int(matches[0])
    else:
        return -1
